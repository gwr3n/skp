package skp.utililities.hash;

import java.nio.charset.StandardCharsets;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;

public class SHA {
   private static String bytesToHex(byte[] hash) {
      StringBuffer hexString = new StringBuffer();
      for (int i = 0; i < hash.length; i++) {
         String hex = Integer.toHexString(0xff & hash[i]);
         if(hex.length() == 1) hexString.append('0');
         hexString.append(hex);
      }
      return hexString.toString();
   }

   public static String generateSHA256(String msg) {
      MessageDigest digest;
      try {
         digest = MessageDigest.getInstance("SHA-256");
         byte[] encodedhash = digest.digest(msg.getBytes(StandardCharsets.UTF_8));
         return bytesToHex(encodedhash);
      } catch (NoSuchAlgorithmException e) {
         // TODO Auto-generated catch block
         e.printStackTrace();
      }
      return null;
   }
}
