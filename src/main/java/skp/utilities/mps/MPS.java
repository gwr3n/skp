package skp.utilities.mps;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.nio.file.StandardCopyOption;

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
   
}
