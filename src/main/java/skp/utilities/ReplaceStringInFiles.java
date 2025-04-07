package skp.utilities;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.stream.Stream;

public class ReplaceStringInFiles {

    public static void main(String[] args) {
        String directoryPath = "./batch/normal";
        String target = "solved_generic_distribution_instances_MILP";
        String replacement = "solved_normal_instances_DCG";

        try {
            traverseAndReplace(directoryPath, target, replacement);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static void traverseAndReplace(String directoryPath, String target, String replacement) throws IOException {
        try (Stream<Path> paths = Files.walk(Paths.get(directoryPath))) {
            paths.filter(Files::isRegularFile).forEach(path -> {
                try {
                    replaceInFile(path, target, replacement);
                } catch (IOException e) {
                    e.printStackTrace();
                }
            });
        }
    }

    public static void replaceInFile(Path path, String target, String replacement) throws IOException {
        File file = path.toFile();
        String content = new String(Files.readAllBytes(file.toPath()));
        content = content.replace(target, replacement);
        Files.write(file.toPath(), content.getBytes(), StandardOpenOption.TRUNCATE_EXISTING);
    }
}
