import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

public class FoldAllCreator {
	
	static final String currentVersion = "java -jar csvscraper1.4.jar";
	
	/*
	 * Either
	 * openacc
	 * cpu2017
	 * cpu2006
	 * mpi2007
	 * omp2012
	 */
	static String testType = "cpu2006";
	
	public static void main(String[] args) throws IOException{
		String s = "";
		
		FileWriter writer = new FileWriter("foldAllFiles.sh");
		
		writer.append(getBeginning());
		writer.append(fileArguments());
		
		writer.close();
		
	}
	public static String getBeginning(){
		return currentVersion + " " + testType + " ";
	}
	public static String fileArguments(){
		String output = "";
		File[] files = new File(System.getProperty("user.dir")+"\\raw").listFiles();
		int i = 0;
		for(File f : files){
			output+=f.getName() + ((i==f.length() -1) ? "" : " ");
			System.out.println(f.getName());
			i++;
		}
		
		
		
		return output;
	}

}
