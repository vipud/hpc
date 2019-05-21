import java.io.FileWriter;
import java.io.IOException;

public class LinkGenerator {
/*
 * 
 * https://www.spec.org/cpu2017/results/res2017q2/cpu2017-20161026-00009.csv
 */
	
	static final String beg = "https://www.spec.org/cpu2017/results/res";
	static final String test = "cpu2017-";
	
	public static void main(String[] args) throws IOException{
		FileWriter writer = new FileWriter("download-links");
		for(int year = 2017;year < 2020;year++){
			for(int q = 1; q < 5;q++){
				String toAdd = "" + beg;
				toAdd += "" + year + "q" + q + "/";
			}
		}
		
		
		
		writer.close();
		
	}
	
}
