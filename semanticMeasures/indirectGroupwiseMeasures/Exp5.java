import java.util.Set;
import org.openrdf.model.URI;
import slib.graph.algo.utils.GAction;
import slib.graph.algo.utils.GActionType;
import slib.graph.algo.utils.GraphActionExecutor;
import slib.graph.io.conf.GDataConf;
import slib.graph.io.loader.GraphLoaderGeneric;
import slib.graph.io.util.GFormat;
import slib.graph.model.graph.G;
//import org.openrdf.model.vocabulary.RDFS;
import slib.graph.model.impl.graph.memory.GraphMemory;
import slib.graph.model.impl.repo.URIFactoryMemory;
import slib.graph.model.repo.URIFactory;
import slib.graph.algo.accessor.GraphAccessor;
import slib.sml.sm.core.engine.SM_Engine;
import slib.sml.sm.core.metrics.ic.utils.IC_Conf_Topo;
import slib.sml.sm.core.metrics.ic.utils.ICconf;
import slib.sml.sm.core.utils.SMConstants;
import slib.sml.sm.core.utils.SMconf;
import slib.utils.ex.SLIB_Exception;

//import exp5.Gene; 
import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.DataOutputStream;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintWriter;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.Random;

import java.io.IOException;
import java.io.InputStreamReader;
import java.io.ObjectOutputStream;
import java.io.OutputStream;
import java.util.Arrays;
import java.util.Iterator;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.logging.Level;
import java.util.logging.Logger;
import slib.utils.ex.SLIB_Ex_Critic;

public class Exp5 {
  public static void main(String[] params) throws SLIB_Exception{
      String gM_indirect_string = "SIM_GROUPWISE_AVERAGE";
	  String pM_string = "SIM_PAIRWISE_DAG_NODE_FEATURE_TVERSKY_RATIO_MODEL";
	  String IC_string = "ICI_PROB_OCCURENCE_PROPAGATED";
      
    //1. CREATE GO ONTOLOGY
      //Access to the in-memory URI Factory.__________________________________
      URIFactory factory = URIFactoryMemory.getSingleton();
      
      //______________________________________________________________________
      URI graph_uri = factory.getURI("http://go/");
      factory.loadNamespacePrefix("GO", graph_uri.toString());
      G graph = new GraphMemory(graph_uri);

      //Load OBO file to graph________________________________________________
      //"gene_ontology_ext.obo"
      GDataConf goConf = new GDataConf(GFormat.OBO, "../gene_ontology_ext.obo");
      GraphLoaderGeneric.populate(goConf, graph);

      //Add virtual root for 3 subontologies__________________________________
      URI virtualRoot = factory.getURI("http://go/virtualRoot");
      graph.addV(virtualRoot);
      GAction rooting = new GAction(GActionType.REROOTING);
      rooting.addParameter("root_uri", virtualRoot.stringValue());
      GraphActionExecutor.applyAction(factory, rooting, graph);

      //Vertices excluding gene products______________________________________
      Set<URI> Vertices = GraphAccessor.getClasses(graph);
     
    //2. ANNOTATE 5000 GENE PRODUCTS FROM ANNOTATION FILE 
      //Constants_____________________________________________________________
      int geneNum = 5000;                 //Number of gene products
      int annotSizesNum = 5;              //Number of different annotation sizes
      int gN_div_aN = geneNum/annotSizesNum;
      Gene[] gene = new Gene[geneNum];               //Allocate array
      for (int i = 0; i < geneNum; i++) {
        gene[i] = new Gene(i + 1, new LinkedHashSet());    //Initialize array of gene products
      }
            
      FileInputStream fis = null;
      BufferedReader reader = null;
      try {
        fis = new FileInputStream("../annotations_plain.txt");
        reader = new BufferedReader(new InputStreamReader(fis));
        String line;  
        
        //i) 1000 gene products with annotation size 1
        for (int i = 0; i < geneNum/annotSizesNum; i++) {
          line = reader.readLine();
          for (URI j : Vertices) {
            if (j.toString().equals(line)) {
              System.out.println("correct annot" + j.toString());
              gene[i].annot.add(j); 
            }
          }                  
        }               
        //ii) 1000 gene products with annotation size 10
        for (int i = gN_div_aN; i < 2*gN_div_aN; i++) {
          for (int j = 0; j < 10; j++) {
            line = reader.readLine();
              for (URI k : Vertices) {
                if (k.toString().equals(line)) {
                  gene[i].annot.add(k); 
                }
              }  
           }
        }  
        //iii) 1000 gene products with annotation size 50
        for (int i = 2*gN_div_aN; i < 3*gN_div_aN; i++) {
          for (int j = 0; j < 50; j++) {   
            line = reader.readLine();
            for (URI k : Vertices) {
              if (k.toString().equals(line)) {
                gene[i].annot.add(k); 
              }
            }  
          }
        } 
        //iv) 1000 gene products with annotation size 100
        for (int i = 3*gN_div_aN; i < 4*gN_div_aN; i++) {
          for (int j = 0; j < 100; j++) {
            line = reader.readLine();
            for (URI k : Vertices) {
              if (k.toString().equals(line)) {
                gene[i].annot.add(k); 
              }
            }  
          }
        } 
        
        //v) 1000 gene products with annotation size 1000
        for (int i = 4*gN_div_aN; i < geneNum; i++) {
          for (int j = 0; j < 1000; j++) {
           line = reader.readLine();
           for (URI k : Vertices) {
              if (k.toString().equals(line)) {
                gene[i].annot.add(k); 
              }
            }  
          } 
        } 
            
      } 
      catch (FileNotFoundException ex) {
        Logger.getLogger(Exp5.class.getName()).log(Level.SEVERE, null, ex);
      } 
      catch (IOException ex) {
        Logger.getLogger(Exp5.class.getName()).log(Level.SEVERE, null, ex);   
      }
      finally {
        try {
          reader.close();
          fis.close();
        }
        catch (IOException ex) {
          Logger.getLogger(Exp5.class.getName()).log(Level.SEVERE, null, ex);
        }
      }
           
    //3. CREATE SETS OF ALL POSSIBLE PAIR AND GROUP MEASURES     
      SM_Engine engine = new SM_Engine(graph);                                  //Engine based on graph  
      SMconf gM_indirect = new SMconf(gM_indirect_string,gM_indirect_string);         //Direct groupwise measures
      SMconf pM = new SMconf(pM_string,pM_string);
      pM.setICconf(new IC_Conf_Topo(IC_string, IC_string));
	  
    //4. CALCULATE ALL MEASURES
      double[] similarities_indirect = new double[geneNum*geneNum];   
      //Direct groupwise similarity in parallel.
      Arrays.parallelSetAll( similarities_indirect, i -> {
        try {
          return engine.compare(gM_indirect, pM, gene[i/geneNum].annot, gene[i%geneNum].annot);
        }
        catch (SLIB_Ex_Critic ex) {
            return -2.0;
          //Logger.getLogger(Exp3.class.getName()).log(Level.SEVERE, null, ex);
        }
      });
         
    //5. NEATLY WRITE GENE PRODUCT SIMILARITIES IN BINRY TO INDIVIDUAL FILES    
      try {
        OutputStream os = new FileOutputStream( "similarities_" + IC_string + gM_indirect_string + pM_string + ".bin" );
        DataOutputStream dos = new DataOutputStream( os );
        for (int i = 0; i < geneNum*geneNum; i++) {
          dos.writeDouble(similarities_indirect[i]);
          //System.out.println(similarities_direct[i]);
        } 
        dos.close();
      } catch (IOException ex) {
        Logger.getLogger(Exp5.class.getName()).log(Level.SEVERE, null, ex);
      }   
        
      
  }   
}

