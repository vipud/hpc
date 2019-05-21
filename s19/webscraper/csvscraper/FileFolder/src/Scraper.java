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
	//openACC
	/*static final String[] keywords = {"\"Hardware Model:\"",
									  "SPECaccel_acc_base",
									  "\"CPU Name\"",
									  "\"CPU(s) enabled\"",
									  "\"Accel Model Name\"",
									  "Memory",
									  "\"Operating System\"",
									  "\"Other Software\"",
									  "Compiler"};*/
	static final String[] keywords_openacc = {"\"Hardware Model:\"",
			  								  "SPECaccel_acc_base",
			  								  "\"CPU Name\"",
			  								  "\"CPU(s) enabled\"",
			  								  "\"Accel Model Name\"",
			  								  "Memory",
			  								  "\"Operating System\"",
			  								  "\"Other Software\"",
	                                          "Compiler"};
	static final String[] keywords_cpu2017 = {"\"Hardware Model:\"",
											  "SPECrate2017_fp_base",
											  "\"CPU Name\"",
											  "Enabled",
											  "Memory",
											  "OS",
											  "Compiler"};
	static final String[] keywords_cpu2006 = {"\"Hardware Model:\"",
					                          "SPECint_base2006",
					                          "\"CPU Name\"",
					                          "\"CPU(s) enabled\"",					                         
					                          "Memory",
					                          "\"Operating System\"",
					                          "\"Other Software\"",
											  "Compiler"};
	static final String[] keywords_omp2012 = {"\"Hardware Model:\"",
											  "SPECompG_base2012",
			                                  "\"CPU Name\"",
			                                  "\"CPU(s) enabled\"",
			                                  "Memory",
			                                  "\"Operating System\"",
	                                          "Compiler"};
	static final String[] keywords_mpi2007 = {"Model",
											  "SPECmpiL_base2007",
											  "\"Chips enabled\"",
											  "\"Cores enabled\"",
											  "\"Cores per chip\"",
											  "\"Threads per core\"",
											  "Memory",
											  "\"Operating System\"",
											  "\"Other Software\"",};
	static String[][] keyWords = {keywords_openacc,keywords_cpu2017,keywords_cpu2006,keywords_omp2012,keywords_mpi2007};
	 
	static String currentTest = "undefined";
	static String[] keywords;
	
	public static void main(String[] args) throws IOException{
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
		int i = 0;
		for(String s : args){
			if(isCSV(s)){
				if(i == 0){
					keywords = keywords_openacc;
				}
				System.out.println("Folding " + s);
				currentTest=s;
				stringToCSV(buildString(csvToString(s)));
			}else{
				i++;
				assignKeywords(s);
			}
		}
	}
	public static void assignKeywords(String dataset){
		/*
		 * openacc
		 * cpu2017
		 * cpu2006
		 * mpi2007
		 * omp2012
		 */
		if(dataset.equals("openacc")){
			keywords = keywords_openacc;
		}else if(dataset.equals("cpu2006")){
			keywords = keywords_cpu2006;
		}else if(dataset.equals("cpu2017")){
			keywords = keywords_cpu2017;
		}else if(dataset.equals("mpi2007")){
			keywords = keywords_mpi2007;
		}else if(dataset.equals("omp2012")){
			keywords = keywords_omp2012;
		}
	}
	public static boolean isCSV(String s){
		if(s.length() >= 4){
			return s.substring(s.length()-4).equals(".csv");
		}
		return false;
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
		return arr[0];
	}

}
