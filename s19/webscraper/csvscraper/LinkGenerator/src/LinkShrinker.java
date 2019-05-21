import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashSet;

public class LinkShrinker {
	
	public static void main(String[] args) throws IOException{
		HashSet<String> set = new HashSet<>();
		BufferedReader br = new BufferedReader(new FileReader(new File("rawlinkscpu2006.txt")));
		FileWriter writer = new FileWriter("links-extra.txt");
		String line = "";
		while((line = br.readLine()) != null){
			if(line.substring(line.length()-4).equals(".csv")){
				if(set.add(line)){
					writer.append(line+System.lineSeparator());
					System.out.println(line);
				}
			}
		}
		writer.close();
	}

}
