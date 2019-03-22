import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;

public class Scraper {
	/*
	 * this is so easy to make but so terribly inefficient 
	 * as far as space complexity is concerned....
	 * but if it works, it works, and I dont think the .csv files are 
	 * that big, and this won't take that long to run on all of them
	 */
	static final String[] keywords = {"SPECaccel_acc_base",
									  "\"CPU Name\"",
									  "\"CPU(s) enabled\"",
									  "\"Accel Model Name\"",
									  "Memory",
									  "\"Operating System\""};
	
	public static void main(String[] args) throws IOException{
		/*String s = "blasdhkasdfh\nSPECaccel_acc_base,1.740416,,,,\n"
				+ "etckjasdhfksdaf,,,dshfas,jasdfkj,\n\n,ashfasdhf,\n"
				+ " \"Accel Model Name\",\"Tesla K20\" ";
		System.out.println(getNext(s,"SPECaccel_acc_base"));
		System.out.println(getNext(s,"\"Accel Model Name\""));
		stringToCSV("does,this,work?");*/
		
		String test = csvToString("src/accel-20140303-00010.csv");
		System.out.println(buildString(test));
	}
	/*
	 * Given a .csv file, return it as a big string.
	 */
	public static String csvToString(String fname) throws IOException{
		String contents = new String(Files.readAllBytes(Paths.get(fname)));
		return contents;
	}
	
	public static String buildString(String wholeFile){
		
		String output = "";
		for(int i = 0; i < keywords.length;i++){
			output+=getNext(wholeFile,keywords[i])+(i==keywords.length-1 ? "" : ",");
		}
		
		return output;
	}
	/*
	 * given our thing we get returned from buildString, output a .csv file somewhere.
	 * boom done.
	 */
	public static void stringToCSV(String s) throws IOException{
		String[] arr = s.split(",(?=(?:[^\"]*\"[^\"]*\")*[^\"]*$)", -1);
		String appender = "";
		if(arr.length >2){
			appender+=arr[1];
		}
		FileWriter writer = new FileWriter(appender + "-output.csv");
		writer.append(s);
		writer.close();
	}
	/*
	 * Given Strings s1 and s2, return a string that is the next
	 * "entry" contained within s1 that follows after the string
	 * s2.
	 * 
	 * return null if s2 is not found within s1
	 */
	public static String getNext(String s1,String s2){
		int x = s1.indexOf(s2) + s2.length();
		if(x==s2.length()-1)
			return null;
		String[] arr = s1.substring(x).split(",");
//		System.out.println(Arrays.toString(arr));
		return arr[1];
	}

}
