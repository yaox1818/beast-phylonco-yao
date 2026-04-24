package phylonco.beast.evolution.readcountmodel;

import beast.base.core.Function;
import beast.base.inference.MCMC;
import beast.base.inference.Operator;
import beast.base.inference.State;
import beast.base.inference.distribution.Dirichlet;
import beast.base.inference.distribution.ParametricDistribution;
import beast.base.inference.distribution.Prior;
import beast.base.inference.parameter.RealParameter;
import beast.base.util.Randomizer;
import beastfx.app.beast.BeastMCMC;
import org.apache.commons.math3.stat.StatUtils;
import org.junit.jupiter.api.Test;
import org.xml.sax.SAXException;

import javax.xml.parsers.ParserConfigurationException;
import java.io.File;
import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;

import static beast.pkgmgmt.BEASTClassLoader.addServices;
import static org.junit.jupiter.api.Assertions.assertEquals;

public class ExcangeGibbsOperatorTest {

    @Test
    public void test() throws Exception {
        Path dir = Path.of("src","test", "resources");
        File xmlFile = Paths.get(dir.toString(),"ExchangeGibbsOperatorTest.xml").toFile();
        final BeastMCMC beastMCMC = new BeastMCMC();
        final List<String> MCMCargs = new ArrayList<>();
        MCMCargs.add("-overwrite");
        MCMCargs.add(xmlFile.getAbsolutePath());
        beastMCMC.parseArgs(MCMCargs.toArray(new String[0]));
        beastMCMC.run();
    }

}
