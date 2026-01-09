package skp.utilities;

import java.io.IOException;
import java.nio.file.*;
import java.util.stream.Stream;

public class RenameFiles {

    public static void main(String[] args) {
        String directoryPath = "./batch/normal";
        String target = "solved_generic_distribution_instances_MILP";
        String replacement = "solved_normal_instances_DCG";

        try {
            traverseAndRename(directoryPath, target, replacement);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static void traverseAndRename(String directoryPath, String target, String replacement) throws IOException {
        try (Stream<Path> paths = Files.walk(Paths.get(directoryPath))) {
            paths.filter(Files::isRegularFile).forEach(path -> {
                try {
                    renameFile(path, target, replacement);
                } catch (IOException e) {
                    e.printStackTrace();
                }
            });
        }
    }

    public static void renameFile(Path path, String target, String replacement) throws IOException {
        String fileName = path.getFileName().toString();
        if (fileName.contains(target)) {
            String newFileName = fileName.replace(target, replacement);
            Files.move(path, path.resolveSibling(newFileName));
        }
    }
}
