package skp.utililities.gson;

import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import com.google.gson.JsonIOException;
import com.google.gson.stream.JsonReader;

public class GSONUtility<T> {
   
   public static <T> String printInstanceAsJSON(T instance){
      Gson gson = new GsonBuilder().serializeSpecialFloatingPointValues().setPrettyPrinting().create();
      return gson.toJson(instance);
   }
   
   public static <T> void saveInstanceToJSON(T instance, String fileName){
      Gson gson = new GsonBuilder().serializeSpecialFloatingPointValues().setPrettyPrinting().create();
      try {
         FileWriter fw = new FileWriter(fileName);
         gson.toJson(instance, fw);
         fw.close();
      } catch (JsonIOException e) {
         // TODO Auto-generated catch block
         e.printStackTrace();
      } catch (IOException e) {
         // TODO Auto-generated catch block
         e.printStackTrace();
      }
   }
   
   public static <T> T retrieveJSONInstance(String fileName, Class<T> type){
      Gson gson = new Gson();
      JsonReader reader;
      try {
         FileReader fr = new FileReader(fileName);
         reader = new JsonReader(fr);
         T instance = gson.fromJson(reader, type);
         fr.close();
         return instance;
      } catch (IOException e) {
         e.printStackTrace();
         return null;
      }
   }
}
