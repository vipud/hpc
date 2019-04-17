import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Arrays;

public class Scraper {
	/*
	 * this is so easy to make but so terribly inefficient 
	 * as far as space complexity is concerned....
	 * but if it works, it works, and I dont think the .csv files are 
	 * that big, and this won't take that long to run on all of them
	 */
	static final String[] keywords = {"\"Hardware Model:\"",
									  "SPECaccel_acc_base",
									  "\"CPU Name\"",
									  "\"CPU(s) enabled\"",
									  "\"Accel Model Name\"",
									  "Memory",
									  "\"Operating System\""};
	static String currentTest = "undefined";
	
	public static void main(String[] args) throws IOException{
		System.out.println("The JAR works.");
		/*String s = "blasdhkasdfh\nSPECaccel_acc_base,1.740416,,,,\n"
				+ "etckjasdhfksdaf,,,dshfas,jasdfkj,\n\n,ashfasdhf,\n"
				+ " \"Accel Model Name\",\"Tesla K20\" ";
		System.out.println(getNext(s,"SPECaccel_acc_base"));
		System.out.println(getNext(s,"\"Accel Model Name\""));
		stringToCSV("does,this,work?");*/
		
		/*
		currentTest="accel-20140228-00005.csv";
		String test = csvToString("src/accel-20140228-00005.csv");
		System.out.println(buildString(test));
		stringToCSV(buildString(test));
		System.out.println("Done");
		*/
		
		for(String s : args){
			System.out.println("Folding " + s);
			currentTest=s;
			stringToCSV(buildString(csvToString(s)));
		}
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
		output=output.replaceAll("\n", "");
		return output;
	}
	/*
	 * given our thing we get returned from buildString, output a .csv file somewhere.
	 * boom done.
	 */
	public static void stringToCSV(String s) throws IOException{
		System.out.println("We have read the file and folded it into the following string?.");
		s=s.replaceAll("\\R", "");
		System.out.println(s);
		System.out.println("RAHH");
		String[] arr = s.split(",(?=(?:[^\"]*\"[^\"]*\")*[^\"]*$)", -1);
		System.out.println("");
		System.out.println(Arrays.toString(arr));
		
		String appender = currentTest;
		
		s+=",\"" + currentTest.substring(currentTest.indexOf("accel-")+6, currentTest.lastIndexOf(".csv")) + "\"";
		System.out.println("final s:\n" + s + "\nrahhrhashdf");

/*
		if(arr.length >2){
			appender+=arr[1];
		}
		appender = appender.replaceAll("\"", "");
		*/
		FileWriter writer = new FileWriter("scraped-"+appender);
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
		
		String[] arr = s1.substring(x).split(",(?=(?:[^\"]*\"[^\"]*\")*[^\"]*$)", -1)[1].split("\n");
		System.out.println("looking for " + s2 + " and we get the following split");
		System.out.println(Arrays.toString(arr));
		System.out.println("first element is " + arr[0]);
		return arr[0];
	}

}
