import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.stream.Collectors;

import org.jsoup.Jsoup;
import org.jsoup.nodes.Document;
import org.jsoup.nodes.Element;
import org.jsoup.select.Elements;

/**
 * Example program to list links from a URL.
 */

public class ListLinks {


	static ArrayList<String> broken;
    public static void main(String[] args) throws IOException {
       // Validate.isTrue(args.length == 1, "usage: supply url to fetch");
    	//String url = args[0];
    	broken=new ArrayList<String>();
        String url="https://www.spec.org/mpi2007/results/";
  
        HashSet<String> folders = (HashSet<String>) getLinks(url).stream().filter(s -> s.contains("/results/res")).collect(Collectors.toSet());
        HashSet<String> links = new HashSet<>();
        
        for(String s : folders){
        	print(s);
        }
        /*
         * Some web pages don't work? This could be a jsoup version issue.
         * https://www.spec.org/cpu2006/results/res2014q4/
		   https://www.spec.org/cpu2006/results/res2012q3/
		   https://www.spec.org/cpu2006/results/res2014q1/
		   https://www.spec.org/cpu2006/results/res2014q3/
		   https://www.spec.org/cpu2006/results/res2011q2/
		   https://www.spec.org/cpu2006/results/res2010q3/
		   https://www.spec.org/cpu2006/results/res2017q2/
		   https://www.spec.org/cpu2006/results/res2016q3/
		   https://www.spec.org/cpu2006/results/res2015q4/
		   https://www.spec.org/cpu2006/results/res2007q4/
		   https://www.spec.org/cpu2006/results/res2009q4/
		   https://www.spec.org/cpu2006/results/res2012q2/
		   https://www.spec.org/cpu2006/results/res2014q2/
		   https://www.spec.org/cpu2006/results/res2012q4/
	       https://www.spec.org/cpu2006/results/res2011q1/
         */
     
        for(String s : folders){
        	HashSet<String> toAdd;
        	if((toAdd = getLinks(s))!=null)
        		links.addAll(toAdd);
        }
        links = (HashSet<String>) links.stream().filter(s -> s.substring(s.length()-4).equals(".csv")).collect(Collectors.toSet());
        
        
        for(String s : links){
        	System.out.println(s);
        }
        sendToFile(links);
        System.out.println("The following links are broken in some way");
        for(String s : broken){
        	System.out.println(s);
        }
        
    }
    public static void sendToFile(HashSet<String> strings) throws IOException{
    	FileWriter fw = new FileWriter("links-cpu2006.txt");
    	for(String s : strings){
    		fw.append(s + System.lineSeparator());
    	}
    	fw.close();
    }
    /*
     * Given a url, return a HashSet of all Links.
     */
    public static HashSet<String> getLinks(String url) throws IOException{
    	HashSet<String> output = new HashSet<>();
    	print("Fetching %s...", url);
    	Document doc;
    	try{
    		doc = Jsoup.connect(url).get();
    	}catch(IOException e){
    		broken.add(url);
    		return null;
    	}
    	Elements links = doc.select("a[href]");
    	for(Element link : links){
    		output.add(link.attr("abs:href"));
    	}
    	return output;
    }

    private static void print(String msg, Object... args) {
        System.out.println(String.format(msg, args));
    }

}