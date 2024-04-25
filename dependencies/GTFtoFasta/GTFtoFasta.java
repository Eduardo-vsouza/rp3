import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.zip.GZIPOutputStream;

/**
 * GTFtoFasta takes a GTF file and a genomic sequence file that it corresponds to and outputs:
 * 1) two versions of the nucleotide sequences corresponding to GTF transcripts (unique transcript id)
 * 2) the 3frame translated peptide sequences corresponding to 1)
 * 3) a split up version of the GTF file corresponding to the peptides in 2)
 * @author Maxim Shokhirev. 
 * Copyright (c)  2016-18   The Salk Institute for Biological Studies.
	All Rights Reserved
	 
	Permission to copy, modify and distribute any part of this GTFtoFasta for educational, research and non-profit purposes, without fee, and without a written agreement is hereby granted, provided that the above copyright notice, this paragraph and the following three paragraphs appear in all copies.
	Those desiring to incorporate this GTFtoFasta into commercial products or use for commercial purposes should contact the Technology Transfer Office, Salk Institute for Biological Studies, La Jolla, 10010 N Torrey Pines Rd., La Jolla, CA 92037, Ph: (858) 453-4100.
	IN NO EVENT SHALL THE SALK INSTITUTE FOR BIOLOGICAL STUDIES BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS, ARISING OUT OF THE USE OF THIS GTFtoFASTA, EVEN IF THE SALK INSTITUTE FOR BIOLOGICAL STUDIES HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
	THE GTFtoFASTA PROVIDED HEREIN IS ON AN "AS IS" BASIS, AND THE SALK INSTITUTE FOR BIOLOGICAL STUDIES HAS NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.  THE SALK INSTITUTE FOR BIOLOGICAL STUDIES MAKES NO REPRESENTATIONS AND EXTENDS NO WARRANTIES OF ANY KIND, EITHER IMPLIED OR EXPRESS, 
	INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE, OR THAT THE USE OF THE GTFtoFASTA WILL NOT INFRINGE ANY PATENT, TRADEMARK OR OTHER RIGHTS.
 */

public class GTFtoFasta {

	HashMap<String,HashMap<String,LinkedList<Integer>>> transcripts = new HashMap<String,HashMap<String,LinkedList<Integer>>>();
	HashMap<String,char[]> chrs = new HashMap<String,char[]>();

	public GTFtoFasta(String[] args) throws IOException {
		if(args.length < 2)
		{
			System.err.println("Usage: GTFtoFasta transcripts.gtf genome.fa [MsOnly]\n"
					+ "Uses genome sequence to produce 3-frame translated transcripts split up by stop codons\n"
					+ "Outputs the transcript dna sequences to *.nuc.gz\n"
					+ "Outputs the transcript peptide sequence to *.pep.gz\n"
					+ "Outputs the updated gtf file to *_ORFs.gtf.gz"
					+ "\nif the 3rd arg is present, will output Met to Stop when possible instead of Stop-Stop");
			System.exit(1);
		}
		if(args.length > 2)
			System.err.println("Will output Met-Stop where possible instead of Stop-Stop");
		SuperScanner gtfScanner = new SuperScanner(args[0]);
		while(gtfScanner.hasMore())
		{
			String line = gtfScanner.getLine();
			String[] split = line.split("\t");
			if(split.length < 3) continue;
			String chr = split[0];
			if(split[2].compareTo("exon")==0)
			{
				if(transcripts.get(chr) == null) transcripts.put(chr,new HashMap<String,LinkedList<Integer>>());
				String transcript = new String(split[6]+line.substring(line.indexOf("transcript_id")+15));
				transcript = transcript.substring(0,transcript.indexOf("\""));
				if(transcripts.get(chr).get(transcript) == null)
					transcripts.get(chr).put(transcript, new LinkedList<Integer>());
				transcripts.get(chr).get(transcript).add(Integer.parseInt(split[3]));
				transcripts.get(chr).get(transcript).add(Integer.parseInt(split[4]));
				
			}
		}
		for(String chr: transcripts.keySet())
		{
			for(String id: transcripts.get(chr).keySet())
			{
				Collections.sort(transcripts.get(chr).get(id));
			}
		}
		//now all transcritps are sorted by their exon order from left to right along the genome
		gtfScanner.close();
		//lets sort all of the transcripts such that their exons are in order
		
		//now that we have read in all of the transcripts lets load in the entire genome
		SuperScanner fa = new SuperScanner(args[1]);
		StringBuilder sb = new StringBuilder(30000000);
		String chr = null;
		while(fa.hasMore())
		{
			String line = fa.getLine();
			if(line.startsWith(">"))
			{
				if(chr != null)
				{
					chrs.put(chr,sb.toString().toCharArray());
					System.err.println("Loaded "+chr);
				}
				chr = line.substring(1);
				sb = new StringBuilder(30000000);
			}
			else
			{
				sb.append(line);
			}
		}
		if(chr != null)
		{
			chrs.put(chr,sb.toString().toCharArray());
			System.err.println("Loaded "+chr);
			sb = null;
		}
		PrintWriter nucWriter = new PrintWriter(new GZIPOutputStream(new FileOutputStream(args[0].substring(0, args[0].indexOf(".gtf"))+".nuc.gz")));
		PrintWriter splitNucWriter = new PrintWriter(new GZIPOutputStream(new FileOutputStream(args[0].substring(0, args[0].indexOf(".gtf"))+".split_nuc.gz")));
		PrintWriter pepWriter = new PrintWriter(new GZIPOutputStream(new FileOutputStream(args[0].substring(0, args[0].indexOf(".gtf"))+".pep.gz")));
		PrintWriter gtfWriter = new PrintWriter(new GZIPOutputStream(new FileOutputStream(args[0].substring(0, args[0].indexOf(".gtf"))+"_ORFs.gtf.gz")));
		StringBuilder nuc = new StringBuilder(1000000);
		StringBuilder snuc = new StringBuilder(1000000);
		StringBuilder three = new StringBuilder(1000000);
		StringBuilder gtf = new StringBuilder(1000000);
		//now lets go through and for each transcript, stitch together its exon sequences
		int count = 0;
		for(String c: transcripts.keySet())
		{
			count++;
			HashMap<String,LinkedList<Integer>> list = transcripts.get(c);
			for(String t: list.keySet())
			{
				LinkedList<Integer> es = list.get(t);
				LinkedList<Integer> copy = new LinkedList<Integer>();
				copy.addAll(es);
				String seq = "";
				boolean plus = true;
				if(t.startsWith("-"))
					plus = false;
				t=t.substring(1);
				int tstart = -1;
				int tend =0;
				while(es.size() > 0)
				{
					int start = es.removeFirst();
					int end = es.removeFirst();
					if(tstart == -1) tstart = start;
					tend=end;
					seq+=getSequence(c,start-1,end);
				}
				if(!plus)
					seq = revComplement(seq);
				//ok now we want to output the transcript in all 3 frames
				if(seq.length() >= 18)
				{

					for(int frame = 0; frame <=2; frame ++)
					{
						nuc.append(">");
						int actualFrame = (tstart+frame)%3;
						if(!plus) actualFrame = (tend-frame)%3;
						String peptide = seqToAA(seq.substring(frame));
						String[] frags = peptide.split("X");
						String[] nucfrags = new String[frags.length];
						String temp = new String(seq.substring(frame));
						for(int fragi = 0; fragi < frags.length; fragi++)
						{
							nucfrags[fragi]=temp.substring(0,frags[fragi].length()*3);
							if(fragi < frags.length-1)
								temp = temp.substring((frags[fragi].length()+1)*3);
						}
						int fragCount = 0;
						int p = frame;  // pointer in the sequence thus far (use regular frame because that is how sequence offset is determined)
						for(int fragi = 0; fragi < frags.length; fragi++)
						{					
							String frag = frags[fragi];
							String nucFrag = nucfrags[fragi];
							if(args.length > 2 && frag.contains("M"))
							{
								int mindex =frag.indexOf("M");
								p+=mindex*3;
								frag = frag.substring(mindex);
								nucFrag = nucFrag.substring(mindex*3);
							}
							if(frag.length() >= 6) //minimum size
							{
								String annotation = "";
								if(frag.contains("M"))
									annotation="_M";
								int p1 = p;
								int p2 = p+frag.length()*3;
								LinkedList<Integer> points = new LinkedList<Integer>();
								if(plus)
								{
									for(int i = 0; i < copy.size(); i+=2)
									{
										int size = copy.get(i+1)-copy.get(i)+1;
										boolean add_right = true;
										boolean add_left = true;
										if(p1 <= size && p1 >= 0)  //the start is in here somewhere
										{
											int start = copy.get(i)+p1;
											points.add(start);
											add_left = false;  //we are passed this
											p1=-1;
										}
										else
										{
											p1-=size;
											//do nothing since we haven't found the start yet
										}
										if(p2 <= size && p1 < 0)
										{
											//end is in here somewhere
											int end = copy.get(i)+p2;
											if(points.size()%2 == 0)points.add(copy.get(i));
											points.add(end);
											break;  //we are done scanning
										}
										else
										{
											p2-=size;
										}
										if(p1< 0)
										{
											if(add_left)
												points.add(copy.get(i));
											if(add_right)
												points.add(copy.get(i+1));
										}

									}
									//now lets output all of the points
									if(points.size() > 1)
									{
										Collections.sort(points);
										int left_most = points.getFirst();
										int right_most = points.getLast()-1;
										String id = t+"+"+c+":"+left_most+"-"+right_most+"_F:"+actualFrame+"_P:"+fragCount+annotation;
										gtf.append(c+"\tGTF2FASTA\ttranscript\t"+left_most+"\t"+right_most+"\t1000\t+\t.\tgene_id \""+id+"\"; transcript_id \""+id+"\";");
										gtf.append(System.lineSeparator());			
										while(points.size() > 1)
										{
											int left = points.removeFirst();
											int right = points.removeFirst();
											if(points.size() == 0)
												right--;
											gtf.append(c+"\tGTF2FASTA\texon\t"+left+"\t"+right+"\t1000\t+\t.\tgene_id \""+id+"\"; transcript_id \""+id+"\";");
											gtf.append(System.lineSeparator());			
										}

										three.append(">"+t+"+"+c+":"+left_most+"-"+right_most+"_F:"+actualFrame+"_P:"+fragCount+annotation);
										three.append(System.lineSeparator());
										three.append(frag);
										three.append(System.lineSeparator());
										snuc.append(">"+t+"+"+c+":"+left_most+"-"+right_most+"_F:"+actualFrame+"_P:"+fragCount+annotation);
										snuc.append(System.lineSeparator());
										snuc.append(nucFrag);
										snuc.append(System.lineSeparator());
										nuc.append(t+"+"+c+":"+left_most+"-"+right_most+"_F:"+actualFrame+"_P:"+fragCount+annotation+" ");

									}
									p+=(frag.length()+1)*3; // to account for the stop codon

								}
								else  //count from the end
								{									
									for(int i = copy.size()-1; i >0 ; i-=2)  //exons are arranged left to right so we want to scan from right to left here
									{
										int size = copy.get(i)-copy.get(i-1)+1;
										boolean add_right = true;
										boolean add_left = true;
										if(p1 <= size && p1 >= 0)  //the start is in here somewhere
										{
											int start = copy.get(i)-p1;
											points.add(start);
											add_right = false;  //we are passed this
											p1=-1;
										}
										else
										{
											p1-=size;
											//do nothing since we haven't found the start yet
										}
										if(p2 <= size && p1 < 0)
										{
											//end is in here somewhere
											int end = copy.get(i)-p2;
											if(points.size()%2 == 0) points.add(copy.get(i));
											points.add(end);
											break;  //we are done scanning
										}
										else
										{
											p2-=size;
										}
										if(p1 < 0)
										{
											if(add_left)
												points.add(copy.get(i-1));
											if(add_right)
												points.add(copy.get(i));
										}
									}

									//now lets output all of the points
									if(points.size() > 1)
									{
										Collections.sort(points);
										points.set(0,points.getFirst()+1);
										int left_most = points.getFirst();
										int right_most = points.getLast();
										String id = t+"-"+c+":"+left_most+"-"+right_most+"_F:"+actualFrame+"_P:"+fragCount+annotation;
										gtf.append(c+"\tGTF2FASTA\ttranscript\t"+left_most+"\t"+right_most+"\t1000\t-\t.\tgene_id \""+id+"\"; transcript_id \""+id+"\";");
										gtf.append(System.lineSeparator());			
				
										while(points.size() > 1)
										{
											int left = points.removeFirst();
											int right = points.removeFirst();
											gtf.append(c+"\tGTF2FASTA\texon\t"+left+"\t"+right+"\t1000\t-\t.\tgene_id \""+id+"\"; transcript_id \""+id+"\";");
											gtf.append(System.lineSeparator());
										}
										three.append(">"+t+"-"+c+":"+left_most+"-"+right_most+"_F:"+actualFrame+"_P:"+fragCount+annotation);
										three.append(System.lineSeparator());
										three.append(frag);
										three.append(System.lineSeparator());
										snuc.append(">"+t+"-"+c+":"+left_most+"-"+right_most+"_F:"+actualFrame+"_P:"+fragCount+annotation);
										snuc.append(System.lineSeparator());
										snuc.append(nucFrag);
										snuc.append(System.lineSeparator());
										nuc.append(t+"-"+c+":"+left_most+"-"+right_most+"_F:"+actualFrame+"_P:"+fragCount+annotation+" ");
									}
									p+=(frag.length()+1)*3; // to account for the stop codon
								}
							}
							else //frag < 6
							{
								p+=(frag.length()+1)*3;
							}
							fragCount++;
						}
						nuc.append(System.lineSeparator());
						nuc.append(seq);
						nuc.append(System.lineSeparator());
					}
				}
				if(gtf.length() > 990000)
				{
					gtfWriter.print(gtf.toString());
					gtfWriter.flush();
					gtf = new StringBuilder(1000000);
					nucWriter.println(nuc.toString());
					nucWriter.flush();
					nuc = new StringBuilder(1000000);
					splitNucWriter.println(snuc.toString());
					splitNucWriter.flush();
					snuc = new StringBuilder(1000000);
					pepWriter.print(three.toString());
					pepWriter.flush();
					three = new StringBuilder(1000000);
				}
			}
			System.err.println("Processed "+c+" ("+count+"/"+transcripts.size()+")");
		}
		if(gtf.length() > 0)
		{
			gtfWriter.print(gtf.toString());
			gtfWriter.flush();
		}
		if(three.length() > 0)
		{
			pepWriter.print(three.toString());
			pepWriter.flush();
		}
		if(nuc.length() > 0)
		{
			nucWriter.println(nuc.toString());
			nucWriter.flush();
		}
		if(snuc.length() > 0)
		{
			splitNucWriter.println(snuc.toString());
			splitNucWriter.flush();
		}
		gtfWriter.close();
		nucWriter.close();
		pepWriter.close();
		splitNucWriter.close();
	}
	private String getSequence(String chr, int start,
			int stop) {
		char[] chars = new char[(int) (stop-start)];
		if(chrs.get(chr) != null)
		{
			char[] bytes = chrs.get(chr);
			try{
				for(int i = start; i < stop; i++)
					chars[i-start]=bytes[i];
			}catch(Exception e)
			{ return "!";}
			return new String(chars).toUpperCase();
		}
		else
			return "!";
	}	
	
	private static String revComplement(String seq)
	{
		char[] arr = seq.toCharArray();
		char[] ret = new char[arr.length];
		for(int i = arr.length-1; i >= 0; i--)
 		{
			if(arr[i]=='A')ret[ret.length-i-1]='T';
			else if(arr[i]=='T')ret[ret.length-i-1]='A';
			else if(arr[i]=='C')ret[ret.length-i-1]='G';
			else if(arr[i]=='G')ret[ret.length-i-1]='C';		
		}
		return new String(ret);
	}
	
	HashMap<String,String> nucToPep = buildNucToPep(); 
	private String seqToAA(String seq)
	{		
		char[] arr = seq.toCharArray();
		char[] res = new char[arr.length/3];
		for(int i = 0; i <= seq.length()-3; i+=3 )
		{
			String trip = ""+arr[i]+arr[i+1]+arr[i+2];
			char aa = '#';
			if(nucToPep.get(trip) != null)
				aa = nucToPep.get(trip).charAt(0);
			res[i/3]=aa;
		}
		return new String(res);
	}
	
	private HashMap<String, String> buildNucToPep() {
		// TODO Auto-generated method stub
		HashMap<String,String> map = new HashMap<String,String>();
		map.put("TTT","F");map.put("TCT","S");map.put("TAT","Y");map.put("TGT","C");
		map.put("TTC","F");map.put("TCC","S");map.put("TAC","Y");map.put("TGC","C");
		map.put("TTA","L");map.put("TCA","S");map.put("TAA","X");map.put("TGA","X");
		map.put("TTG","L");map.put("TCG","S");map.put("TAG","X");map.put("TGG","W");
		map.put("CTT","L");map.put("CCT","P");map.put("CAT","H");map.put("CGT","R");
		map.put("CTC","L");map.put("CCC","P");map.put("CAC","H");map.put("CGC","R");
		map.put("CTA","L");map.put("CCA","P");map.put("CAA","Q");map.put("CGA","R");
		map.put("CTG","L");map.put("CCG","P");map.put("CAG","Q");map.put("CGG","R");
		map.put("ATT","I");map.put("ACT","T");map.put("AAT","N");map.put("AGT","S");
		map.put("ATC","I");map.put("ACC","T");map.put("AAC","N");map.put("AGC","S");
		map.put("ATA","I");map.put("ACA","T");map.put("AAA","K");map.put("AGA","R");
		map.put("ATG","M");map.put("ACG","T");map.put("AAG","K");map.put("AGG","R");
		map.put("GTT","V");map.put("GCT","A");map.put("GAT","D");map.put("GGT","G");
		map.put("GTC","V");map.put("GCC","A");map.put("GAC","D");map.put("GGC","G");
		map.put("GTA","V");map.put("GCA","A");map.put("GAA","E");map.put("GGA","G");
		map.put("GTG","V");map.put("GCG","A");map.put("GAG","E");map.put("GGG","G");
		HashMap<String,String> newMap = new HashMap<String,String>();
		for(String key: map.keySet())
		{
			String copy = new String(key);
			String TtoU=copy.replace('T', 'U');
			newMap.put(TtoU,map.get(key));
			newMap.put(key,map.get(key));
			newMap.put(key.toLowerCase(),map.get(key).toLowerCase());
			newMap.put(TtoU.toLowerCase(), map.get(key).toLowerCase());
		}
		return newMap;	
	}


	public static void main(String[] args) {
		// TODO Auto-generated method stub
		try {
			new GTFtoFasta(args);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

}
