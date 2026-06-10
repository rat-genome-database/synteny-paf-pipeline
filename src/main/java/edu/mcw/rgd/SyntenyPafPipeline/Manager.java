// TODO: investigate whether the PAF coordinate columns (query/target start & end) must be
//       0-based half-open or 1-based. The base of synteny_ucsc start_pos1/stop_pos1 is unconfirmed;
//       if those values are 1-based, the start columns need to be emitted as (start - 1) to conform.
package edu.mcw.rgd.SyntenyPafPipeline;

import edu.mcw.rgd.dao.impl.MapDAO;
import edu.mcw.rgd.dao.impl.SyntenyDAO;
import edu.mcw.rgd.datamodel.SpeciesType;
import edu.mcw.rgd.datamodel.SyntenicRegion;
import edu.mcw.rgd.process.Utils;
import edu.mcw.rgd.process.mapping.MapManager;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.springframework.beans.factory.support.DefaultListableBeanFactory;
import org.springframework.beans.factory.xml.XmlBeanDefinitionReader;
import org.springframework.core.io.FileSystemResource;

import java.io.File;
import java.io.FileWriter;
import java.util.*;


/**
 * Created by cdursun on 2/3/2017.
 */
public class Manager {

    private String version;
    Logger logger = LogManager.getLogger("status");

    /**
     * load spring configuration from properties/AppConfigure.xml file
     * and run the pipeline
     * @param args
     * @throws Exception
     */
    public static void main(String[] args) throws Exception {

        LinkedHashMap<String, Integer> assemblies = new LinkedHashMap<String, Integer>();
        assemblies.put("GRCr8",380);
        assemblies.put("mRatBN7.2",372);
        assemblies.put("Rnor_6.0",360);
        assemblies.put("Rnor_5.0",70);
        assemblies.put("RGSC_v3.4",60);
        assemblies.put("GRCh38.p14",38);
        assemblies.put("GRCh37.p13",17);
        assemblies.put("NCBI36",13);
        assemblies.put("GRCm39",239);
        assemblies.put("GRCm38",35);
        assemblies.put("MGSCv37",18);
        assemblies.put("CanFam3.1",631);
        assemblies.put("UU_Cfam_GSD_1.0",637);
        assemblies.put("Dog10K_Boxer_Tasha",633);
        assemblies.put("ROS_Cfam_1.0",634);
        assemblies.put("Sscrofa11.1",911);
        assemblies.put("Sscrofa10.2",910);
        assemblies.put("Chlorocebus_sabeus1.1",1311);
        assemblies.put("Vero_WHO_p1.0",1313);
        assemblies.put("Mhudiblu_PPA_v0",513);
        assemblies.put("panpan1.1",511);
        assemblies.put("HetGla_female_1.0",1410);
        assemblies.put("ChiLan1.0",44);
        assemblies.put("SpeTri2.0",720);


        DefaultListableBeanFactory bf = new DefaultListableBeanFactory();
        new XmlBeanDefinitionReader(bf).loadBeanDefinitions(new FileSystemResource("properties/AppConfigure.xml"));
        Manager manager = (Manager) bf.getBean("main");

        manager.logger.info(manager.getVersion());


        if( args == null || args.length < 1 ){
            System.out.println("");
            System.out.println("            Missing parameter!                  ");
            System.out.println("----------- Run with output directory specified -----------");
            System.exit(0);
        }

        String outputDirectory = args[0];

        Date time0 = Calendar.getInstance().getTime();

        HashMap<Integer,Boolean> processed = new HashMap<Integer,Boolean>();

        try {

            for (Map.Entry<String, Integer> entry: assemblies.entrySet()) {
                String key = entry.getKey();
                Integer value = entry.getValue();

                for (Map.Entry<String,Integer> entry2: assemblies.entrySet()) {

                    String key2 = entry2.getKey();
                    Integer value2 = entry2.getValue();

                    //if (!key.equals(key2)) {
                   // if (MapManager.getInstance().getMap(value).getSpeciesTypeKey() != MapManager.getInstance().getMap(value2).getSpeciesTypeKey()) {

                        if (!processed.containsKey(value2)) {
                            manager.run(key, value, key2, value2, outputDirectory);
                        }
                    //}

                }
                processed.put(value,true);

            }


            // manager.run(assembly1, mapKey1, assembly2, mapKey2, outputDirectory);
        } catch(Exception e) {
            Utils.printStackTrace(e, manager.logger);
            throw e;
        }

        manager.logger.info("=== OK === elapsed time " + Utils.formatElapsedTime(time0.getTime(), System.currentTimeMillis()));
        manager.logger.info("");
    }

    /*
    Col	Type	Description
1	string	Query sequence name
2	int	Query sequence length
3	int	Query start (0-based; BED-like; closed)
4	int	Query end (0-based; BED-like; open)
5	char	Relative strand: "+" or "-"
6	string	Target sequence name
7	int	Target sequence length
8	int	Target start on original strand (0-based)
9	int	Target end on original strand (0-based)
10	int	Number of residue matches
11	int	Alignment block length
12	int	Mapping quality (0-255; 255 for missing)

     */

    public void run(String assembly1, int mapKey1, String assembly2, int mapKey2, String outputDirectory) throws Exception {

        MapDAO mdao = new MapDAO();

        Map<String,Integer> chrLen1 = mdao.getChromosomeSizes(mapKey1);
        Map<String,Integer> chrLen2 = mdao.getChromosomeSizes(mapKey2);

        SyntenyDAO sdao = new SyntenyDAO();

        // getBlocks is scoped to a single backbone chromosome/range, so aggregate over every chromosome of mapKey1
        List<SyntenicRegion> regions = new ArrayList<>();
        for( String chr: chrLen1.keySet() ) {
            regions.addAll(sdao.getBlocks(mapKey1, chr, 0, chrLen1.get(chr), mapKey2));
        }

        //nothing here
        if (regions.size() == 0) {
            return;
        }

        String track = assembly1 + " (" + SpeciesType.getCommonName(MapManager.getInstance().getMap(mapKey1).getSpeciesTypeKey()) + ") - " + assembly2 + " (" +  SpeciesType.getCommonName(MapManager.getInstance().getMap(mapKey2).getSpeciesTypeKey()) + ")";

        FileWriter fw = new FileWriter(new File(outputDirectory + "/" + track +".paf"));

        for (SyntenicRegion sr: regions) {

            Integer queryLen = chrLen1.get(sr.getBackboneChromosome());
            Integer targetLen = chrLen2.get(sr.getChromosome());

            // every PAF length field must be a valid integer; skip regions whose chromosome size is unknown
            if (queryLen == null || targetLen == null) {
                logger.warn("Skipping region with unknown chromosome length [" + track + "]"
                        + " query chr=" + sr.getBackboneChromosome() + " target chr=" + sr.getChromosome());
                continue;
            }

            long residueMatches = sr.getBackboneStop() - sr.getBackboneStart();

            String row = String.join("\t",
                    "Chr" + sr.getBackboneChromosome(),       // 1  query sequence name
                    queryLen.toString(),                      // 2  query sequence length
                    Integer.toString(sr.getBackboneStart()),  // 3  query start
                    Integer.toString(sr.getBackboneStop()),   // 4  query end
                    sr.getOrientation(),                      // 5  relative strand (+/-)
                    "Chr" + sr.getChromosome(),               // 6  target sequence name
                    targetLen.toString(),                     // 7  target sequence length
                    Integer.toString(sr.getStart()),          // 8  target start
                    Integer.toString(sr.getStop()),           // 9  target end
                    Long.toString(residueMatches),            // 10 number of residue matches
                    Long.toString(residueMatches),            // 11 alignment block length
                    "255");                                   // 12 mapping quality

            fw.write(row + "\n");

        }
        fw.close();

    }

    public void setVersion(String version) {
        this.version = version;
    }

    public String getVersion() {
        return version;
    }
}
