package skp.gurobi;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.nio.file.StandardCopyOption;

/* Copyright 2023, Gurobi Optimization, LLC */

/* This example reads a MIP model from a file, solves it and
   prints the objective values from all feasible solutions
   generated while solving the MIP. Then it creates the fixed
   model and solves that model. */

import com.gurobi.gurobi.*;

public class MPS {
   public static void fixOPLerrors(String instanceIdentifier) {
      // CPLEX generates faulty MPS files in which INDICATORS lines start with FF instead of IF
      // The following code amends the faulty MPS files

      // Create a temporary file to store the modified content
      File tempFile = new File("temp.mps");

      try {         
         // Open the original file and the temporary file
         BufferedReader reader = new BufferedReader(new FileReader(instanceIdentifier+".mps"));
         BufferedWriter writer = new BufferedWriter(new FileWriter(tempFile));

         // Read each line from the original file and modify it if needed
         String line = reader.readLine();
         while (line != null) {
            // If the line starts with FF, replace it with IF
            if (line.startsWith(" FF")) {
               line = line.replaceFirst(" FF", " IF");
            }
            // Write the modified line to the temporary file
            writer.write(line + "\n");
            // Read the next line from the original file
            line = reader.readLine();
         }

         // Close the reader and the writer
         reader.close();
         writer.close();

         // Delete the original file and rename the temporary file
         File originalFile = new File(instanceIdentifier+".mps");
         originalFile.delete();
         tempFile.renameTo(originalFile);

         File folder = new File("mps");
         if (!folder.exists()) folder.mkdirs();
         Files.move(Paths.get(instanceIdentifier+".mps"), Paths.get("mps/"+instanceIdentifier+".mps"), StandardCopyOption.ATOMIC_MOVE);  
      } catch (IOException e) {
         // Handle any exceptions
         e.printStackTrace();
      }
   }

   public static void main(String[] args) {

      try {
         GRBEnv env = new GRBEnv();
         ClassLoader classLoader = env.getClass().getClassLoader();
         GRBModel model = new GRBModel(env, classLoader.getResource("./mps/backlog_900_1_5.0_0.1.mps").getFile());
         if (model.get(GRB.IntAttr.IsMIP) == 0) {
            System.out.println("Model is not a MIP");
            System.exit(1);
         }

         model.optimize();

         int optimstatus = model.get(GRB.IntAttr.Status);
         double objval = 0;
         if (optimstatus == GRB.Status.OPTIMAL) {
            objval = model.get(GRB.DoubleAttr.ObjVal);
            System.out.println("Optimal objective: " + objval);
         } else if (optimstatus == GRB.Status.INF_OR_UNBD) {
            System.out.println("Model is infeasible or unbounded");
            return;
         } else if (optimstatus == GRB.Status.INFEASIBLE) {
            System.out.println("Model is infeasible");
            return;
         } else if (optimstatus == GRB.Status.UNBOUNDED) {
            System.out.println("Model is unbounded");
            return;
         } else {
            System.out.println("Optimization was stopped with status = "
                  + optimstatus);
            return;
         }

         /* Iterate over the solutions and compute the objectives */

         System.out.println();
         for (int k = 0; k < model.get(GRB.IntAttr.SolCount); ++k) {
            model.set(GRB.IntParam.SolutionNumber, k);
            double objn = model.get(GRB.DoubleAttr.PoolObjVal);

            System.out.println("Solution " + k + " has objective: " + objn);
         }
         System.out.println();

         /* Create a fixed model, turn off presolve and solve */

         /*GRBModel fixed = model.fixedModel();

      fixed.set(GRB.IntParam.Presolve, 0);

      fixed.optimize();

      int foptimstatus = fixed.get(GRB.IntAttr.Status);

      if (foptimstatus != GRB.Status.OPTIMAL) {
        System.err.println("Error: fixed model isn't optimal");
        return;
      }

      double fobjval = fixed.get(GRB.DoubleAttr.ObjVal);

      if (Math.abs(fobjval - objval) > 1.0e-6 * (1.0 + Math.abs(objval))) {
        System.err.println("Error: objective values are different");
        return;
      }

      GRBVar[] fvars  = fixed.getVars();
      double[] x      = fixed.get(GRB.DoubleAttr.X, fvars);
      String[] vnames = fixed.get(GRB.StringAttr.VarName, fvars);

      for (int j = 0; j < fvars.length; j++) {
        if (x[j] != 0.0) {
          System.out.println(vnames[j] + " " + x[j]);
        }
      }*/

         // Dispose of models and environment
         //fixed.dispose();
         model.dispose();
         env.dispose();

      } catch (GRBException e) {
         System.out.println("Error code: " + e.getErrorCode() + ". "
               + e.getMessage());
      }
   }
}
