package skp.utilities;

import java.io.IOException;
import java.nio.file.*;
import java.util.stream.Stream;

public class DeleteFiles {

    public static void main(String[] args) {
        String directoryPath = "./batch/normal";
        String target = "solved_normal_instances_SAA";

        try {
            traverseAndDelete(directoryPath, target);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static void traverseAndDelete(String directoryPath, String target) throws IOException {
        try (Stream<Path> paths = Files.walk(Paths.get(directoryPath))) {
            paths.filter(Files::isRegularFile).forEach(path -> {
                try {
                    deleteFile(path, target);
                } catch (IOException e) {
                    e.printStackTrace();
                }
            });
        }
    }

    public static void deleteFile(Path path, String target) throws IOException {
        String fileName = path.getFileName().toString();
        if (fileName.contains(target)) {
            Files.delete(path);
        }
    }
}
