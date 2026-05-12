package phylonco.beast.evolution.readcountmodel;

import beast.base.core.Input;
import beast.base.evolution.operator.TreeOperator;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.util.InputUtil;
import beast.base.util.Randomizer;
import mutablealignment.MATreeLikelihood;
import mutablealignment.MutableAlignment;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class WilsonBaldingGibbsOperator extends TreeOperator {

    public Input<MutableAlignment> mutableAlignmentInput = new Input<>(
            "mutableAlignment",
            "mutable alignment containing genotype states",
            Input.Validate.REQUIRED);

    public Input<MATreeLikelihood> maTreeLikelihoodInput = new Input<>(
            "maTreeLikelihood",
            "tree likelihood for mutable alignment",
            Input.Validate.REQUIRED);

    public Input<LikelihoodReadCountModel> likelihoodReadCountModelInput = new Input<>(
            "likelihoodReadCountModel",
            "read count likelihood model",
            Input.Validate.REQUIRED);

    private MutableAlignment mutableAlignment;
    private MATreeLikelihood maTreeLikelihood;
    private LikelihoodReadCountModel likelihoodReadCountModel;

    private int numStates;
    private int numSites;
    private List<int[]> allGSequences;

    @Override
    public void initAndValidate() {
        mutableAlignment = mutableAlignmentInput.get();
        maTreeLikelihood = maTreeLikelihoodInput.get();
        likelihoodReadCountModel = likelihoodReadCountModelInput.get();
        numStates = mutableAlignment.getDataType().getStateCount();
        numSites = mutableAlignment.getSiteCount();

        // Precompute constant sequences: allGSequences.get(g) = [g, g, ..., g]
        allGSequences = new ArrayList<>(numStates);
        for (int g = 0; g < numStates; g++) {
            int[] seq = new int[numSites];
            Arrays.fill(seq, g);
            allGSequences.add(seq);
        }
    }

    @Override
    public double proposal() {
        final Tree tree = (Tree) InputUtil.get(treeInput, this);
        final int internalNodes = tree.getInternalNodeCount();

        if (internalNodes <= 1) {
            return Double.NEGATIVE_INFINITY;
        }

        double oldMinAge, newMinAge, newRange, oldRange, newAge, logHR_wilsonbalding;

        // choose a random leaf node avoiding root
        final int leafNodeCount = tree.getLeafNodeCount();
        Node i = tree.getNode(Randomizer.nextInt(leafNodeCount));
        final Node p = i.getParent();


        // choose another random node to insert i above
        final int nodeCount = tree.getNodeCount();
        Node j;
        Node jP;

        // make sure that the target branch <k, j> is above the subtree being moved
        do {
            j = tree.getNode(Randomizer.nextInt(nodeCount));
            jP = j.getParent();
        } while ((jP != null && jP.getHeight() <= i.getHeight()) || (i.getNr() == j.getNr()));

        // disallow moves that change the root.
        if (j.isRoot() || p.isRoot()) {
            return Double.NEGATIVE_INFINITY;
        }

        assert jP != null;  // j != root tested above
        final int pnr = p.getNr();
        final int jPnr = jP.getNr();
        if ( jPnr == pnr || j.getNr() == pnr || jPnr == i.getNr())
            return Double.NEGATIVE_INFINITY;

        final Node CiP = getOtherChild(p, i);

        final Node PiP = p.getParent();

        newMinAge = Math.max(i.getHeight(), j.getHeight());
        newRange = jP.getHeight() - newMinAge;
        newAge = newMinAge + (Randomizer.nextDouble() * newRange);
        oldMinAge = Math.max(i.getHeight(), CiP.getHeight());
        oldRange = PiP.getHeight() - oldMinAge;
        logHR_wilsonbalding = Math.log(newRange / Math.abs(oldRange));

        if (oldRange == 0 || newRange == 0) {
            // This happens when some branch lengths are zero.
            // If oldRange = 0, hastingsRatio == Double.POSITIVE_INFINITY and
            // node i can be catapulted anywhere in the tree, resulting in
            // very bad trees that are always accepted.
            // For symmetry, newRange = 0 should therefore be ruled out as well
            return Double.NEGATIVE_INFINITY;
        }

        int iIdx = i.getNr();
        // tree leaf order and alignment column order are linked only by taxon
        // name: do NOT assume k is the alignment column index, look it up.
        int iAlignIdx = mutableAlignment.getTaxonIndex(i.getID());

        // ---- Compute logQ_reverse (under current tree T, before exchange) ----
        int[] currentSeq = mutableAlignment.getSiteValuesByTaxon(iAlignIdx);
        double logQ_reverse = computeLogGibbsProb(iIdx, currentSeq);



        // disconnect p
        final Node pP = p.getParent();
        replace(pP, p, CiP);
        // re-attach, first child node to p
        replace(p, CiP, j);
        // then parent node of j to p
        replace(jP, j, p);

        // mark paths to common ancestor as changed
        if( markCladesInput.get() ) {
            Node iup = pP;
            Node jup = p;
            while (iup != jup) {
                if( iup.getHeight() < jup.getHeight() ) {
                    assert !iup.isRoot();
                    iup = iup.getParent();
                    iup.makeDirty(Tree.IS_FILTHY);
                } else {
                    assert !jup.isRoot();
                    jup = jup.getParent();
                    jup.makeDirty(Tree.IS_FILTHY);
                }
            }
        }

        p.setHeight(newAge);

//        iIdx = i.getNr();
//        iAlignIdx = mutableAlignment.getTaxonIndex(i.getID());

        List<Node> pLeaves = collectLeaves(p);
        Node pLeave0 = pLeaves.get(0);
        int pLeave0Idx = pLeave0.getNr();
        int pLeave0AlignIdx = mutableAlignment.getTaxonIndex(pLeave0.getID());
        maTreeLikelihood.getLogProbsForStateSequence(pLeave0Idx, mutableAlignment.getSiteValuesByTaxon(pLeave0AlignIdx));

        List<Node> pPLeaves = collectLeaves(pP);
        Node pPLeave0 = pPLeaves.get(0);
        int pPLeave0Idx = pPLeave0.getNr();
        int pPLeave0AlignIdx = mutableAlignment.getTaxonIndex(pPLeave0.getID());
        maTreeLikelihood.getLogProbsForStateSequence(pPLeave0Idx, mutableAlignment.getSiteValuesByTaxon(pPLeave0AlignIdx));

        // ---- Gibbs resample + compute logQ_forward (under new tree T') ----
        double[][] treeLogProb = new double[numStates][];
        double[][] rcLogLik = new double[numStates][];
        for (int g = 0; g < numStates; g++) {
            // clone: getLogProbsForStateSequence returns a reference to a shared internal buffer
            treeLogProb[g] = maTreeLikelihood.getLogProbsForStateSequence(iIdx, allGSequences.get(g)).clone();
            rcLogLik[g] = likelihoodReadCountModel.sequenceLogLikelihood(iIdx, allGSequences.get(g));
        }

        int[] newSeq = new int[numSites];
        double logQ_forward = 0.0;
        for (int s = 0; s < numSites; s++) {
            double[] logP = new double[numStates];
            for (int g = 0; g < numStates; g++) {
                logP[g] = treeLogProb[g][s] + rcLogLik[g][s];
            }
            double[] probs = normalizeLogProbs(logP);
            newSeq[s] = sampleFromProbabilities(probs);
            logQ_forward += Math.log(probs[newSeq[s]]);
        }
        mutableAlignment.setSiteValuesByTaxon(iAlignIdx, newSeq);

        // ---- Return compound Hastings ratio ----
        return logHR_wilsonbalding + logQ_reverse - logQ_forward;
    }

    /**
     * Compute log q_G(seq | T, A_{-k}) = sum_s log q(seq[s] | T, A_{-k}, s).
     *
     * Uses the current tree state in maTreeLikelihood (call before or after
     * exchange as appropriate).
     */
    private double computeLogGibbsProb(int k, int[] seq) {
        double[][] treeLogProb = new double[numStates][];
        double[][] rcLogLik = new double[numStates][];
        for (int g = 0; g < numStates; g++) {
            // clone: getLogProbsForStateSequence returns a reference to a shared internal buffer
            treeLogProb[g] = maTreeLikelihood.getLogProbsForStateSequence(k, allGSequences.get(g)).clone();
            rcLogLik[g] = likelihoodReadCountModel.sequenceLogLikelihood(k, allGSequences.get(g));
        }

        double logProb = 0.0;
        for (int s = 0; s < numSites; s++) {
            double[] logP = new double[numStates];
            for (int g = 0; g < numStates; g++) {
                logP[g] = treeLogProb[g][s] + rcLogLik[g][s];
            }
            double logZ = logSumExp(logP);
            logProb += logP[seq[s]] - logZ;
        }
        return logProb;
    }

    private double logSumExp(double[] logP) {
        double max = Double.NEGATIVE_INFINITY;
        for (double v : logP) {
            if (v > max) max = v;
        }
        double sum = 0.0;
        for (double v : logP) {
            sum += Math.exp(v - max);
        }
        return max + Math.log(sum);
    }

    private double[] normalizeLogProbs(double[] logP) {
        double max = Double.NEGATIVE_INFINITY;
        for (double v : logP) {
            if (v > max) max = v;
        }
        double[] probs = new double[logP.length];
        double sum = 0.0;
        for (int j = 0; j < logP.length; j++) {
            probs[j] = Math.exp(logP[j] - max);
            sum += probs[j];
        }
        for (int j = 0; j < probs.length; j++) {
            probs[j] /= sum;
        }
        return probs;
    }

    private int sampleFromProbabilities(double[] probs) {
        double rand = Randomizer.nextDouble();
        double cumulative = 0.0;
        for (int j = 0; j < probs.length; j++) {
            cumulative += probs[j];
            if (rand < cumulative) return j;
        }
        return probs.length - 1;
    }

    /**
     * Collect all leaf nodes under the given node (including the node itself
     * if it is a leaf).
     */
    private List<Node> collectLeaves(Node node) {
        List<Node> leaves = new ArrayList<>();
        if (node.isLeaf()) {
            leaves.add(node);
        } else {
            // Node.getAllLeafNodes(List) adds leaf descendants recursively
            node.getAllLeafNodes(leaves);
        }
        return leaves;
    }
}



