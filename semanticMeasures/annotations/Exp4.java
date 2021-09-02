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

//import exp4.Gene; 
import java.io.BufferedOutputStream;
import java.io.DataOutputStream;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintWriter;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.Random;

import java.io.IOException;
import java.io.ObjectOutputStream;
import java.io.OutputStream;
import java.util.Arrays;
import java.util.Iterator;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.logging.Level;
import java.util.logging.Logger;
import slib.utils.ex.SLIB_Ex_Critic;

public class Exp4 {
  public static void main(String[] params) throws SLIB_Exception{
      
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

      System.out.println(Vertices.size());

    //2. ANNOTATE 5000 GENE PRODUCTS RANDOMLY
      //Constants_____________________________________________________________
      int geneNum = 5000;                 //Number of gene products
      int annotSizesNum = 5;              //Number of different annotation sizes
      int gN_div_aN = geneNum/annotSizesNum;
      System.out.println(gN_div_aN);
      
      //Array of gene products________________________________________________ 
      Gene[] gene = new Gene[geneNum];               //Allocate array
      for (int i = 0; i < geneNum; i++) {
        gene[i] = new Gene(i + 1, new LinkedHashSet());    //Initialize array of gene products
      }
 
      //i) 1000 gene products with annotation size 1
      for (int i = 0; i < gN_div_aN; i++) {
        int labeled = 0;
        while (labeled == 0) {
          //System.out.println("Protein " + gene[i].id + " has annotations:");
          int item = new Random().nextInt(Vertices.size());
          int counter = 0;
          for (URI k : Vertices) {
            if (counter == item && gene[i].annot.contains(k) == false && k.toString().equals(virtualRoot.toString()) == false) {
              //System.out.println("\t" + k);
              gene[i].annot.add(k);
              labeled = 1;
            }
            counter++;
          }
        }
      } 
      
      //ii) 1000 gene products with annotation size 10
      for (int i = gN_div_aN; i < 2*gN_div_aN; i++) { 
        //System.out.println("Protein " + gene[i].id + " has annotations:");
        for (int j = 0; j < 10; j++) {
          int labeled = 0;
          while (labeled == 0) {
            int item = new Random().nextInt(Vertices.size());
            int counter = 0;
            for (URI k : Vertices) {
              if (counter == item && gene[i].annot.contains(k) == false && k.toString().equals(virtualRoot.toString()) == false) {
                //System.out.println("\t" + k);
                gene[i].annot.add(k);
                labeled = 1;
              }
              counter++;
            }
          }
        }
      } 

      //iii) 1000 gene products with annotation size 50
      for (int i = 2*gN_div_aN; i < 3*gN_div_aN; i++) {
        //System.out.println("Protein " + gene[i].id + " has annotations:");
        for (int j = 0; j < 50; j++) {
          int labeled = 0;
          while (labeled == 0) {
            int item = new Random().nextInt(Vertices.size());
            int counter = 0;
            for (URI k : Vertices) {
              if (counter == item && gene[i].annot.contains(k) == false && k.toString().equals(virtualRoot.toString()) == false) {
                //System.out.println("\t" + k);
                gene[i].annot.add(k);
                labeled = 1;
              }
              counter++;
            }
          }
        } 
      }
       
      //iv) 1000 gene products with annotation size 100
      for (int i = 3*gN_div_aN; i < 4*gN_div_aN; i++) {
        //System.out.println("Protein " + gene[i].id + " has annotations:");
        for (int j = 0; j < 100; j++) {
          int labeled = 0;
          while (labeled == 0) {
            int item = new Random().nextInt(Vertices.size());
            int counter = 0;
            for (URI k : Vertices) {
              if (counter == item && gene[i].annot.contains(k) == false && k.toString().equals(virtualRoot.toString()) == false) {
                //System.out.println("\t" + k);
                gene[i].annot.add(k);
                labeled = 1;
              } 
              counter++;
            }
          }   
        }  
      }



      
      //v) 1000 gene products with annotation size 1000
      for (int i = 4*gN_div_aN; i < geneNum; i++) {
        //System.out.println("Protein " + gene[i].id + " has annotations:");
        for (int j = 0; j < 1000; j++) {
          int labeled = 0;
          while (labeled == 0) {
            int item = new Random().nextInt(Vertices.size());
            int counter = 0;
            for (URI k : Vertices) {
              if (counter == item && gene[i].annot.contains(k) == false && k.toString().equals(virtualRoot.toString()) == false) {
                gene[i].annot.add(k);
                labeled = 1;
              }
              counter++;
            }
          }
        } 
      }

for (Gene w : gene) {
  System.out.println(w.id + "\t" + w.annot.size());
} 
      
      
      //Write annotations to file_______________________________________________
      File annot1 = new File("../annotations_plain.txt");
      try {
        PrintWriter printWriter = new PrintWriter(annot1);
        for (int i = 0; i < geneNum; i++) {
          for (URI w : gene[i].annot) {
            //System.out.println(w.toString());
            printWriter.println(w.toString());
          }
        }
        printWriter.close();
      }
      catch (FileNotFoundException ex) {
      }
      
      //Write annotations to file_______________________________________________
      File annot = new File("annotations.txt");
      try {
        //"/data/dragon/yousuff/annotations.txt"
        PrintWriter printWriter = new PrintWriter(annot);
        for (int i = 0; i < geneNum; i++) {
          printWriter.println("Annotations for gene product " + gene[i].id + ":");
          //System.out.println("Annotations for gene product " + gene[i].id + ":");
          for (URI w : gene[i].annot) {
            printWriter.println("\t" + w.toString());
            //System.out.println("\t" + w.toString());
          }
          printWriter.println(" ");
        }
        printWriter.close();
      }
     
      catch (FileNotFoundException ex) {
      }
      
  }   
}
