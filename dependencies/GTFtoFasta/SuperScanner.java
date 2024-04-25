/* Copyright (c)  2016   The Salk Institute for Biological Studies.
	All Rights Reserved
	 
	Permission to copy, modify and distribute any part of this MAPS for educational, research and non-profit purposes, without fee, and without a written agreement is hereby granted, provided that the above copyright notice, this paragraph and the following three paragraphs appear in all copies.
	Those desiring to incorporate this MAPS into commercial products or use for commercial purposes should contact the Technology Transfer Office, The Salk Institute for Biological Studies, La Jolla, 10010 N Torrey Pines Rd., La Jolla, CA 92037, Ph: (858) 453-4100.
	IN NO EVENT SHALL THE SALK INSTITUTE FOR BIOLOGICAL STUDIES BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS, ARISING OUT OF THE USE OF THIS MAPS, EVEN IF THE SALK INSTITUTE FOR BIOLOGICAL STUDIES HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
	THE MAPS PROVIDED HEREIN IS ON AN "AS IS" BASIS, AND THE SALK INSTITUTE FOR BIOLOGICAL STUDIES HAS NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.  THE SALK INSTITUTE FOR BIOLOGICAL STUDIES MAKES NO REPRESENTATIONS AND EXTENDS NO WARRANTIES OF ANY KIND, EITHER IMPLIED OR EXPRESS, 
	INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE, OR THAT THE USE OF THE MAPS WILL NOT INFRINGE ANY PATENT, TRADEMARK OR OTHER RIGHTS.
*/
import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.zip.GZIPInputStream;


public class SuperScanner {

	BufferedReader s = null;
	
	public SuperScanner(String file)
	{
		try{
			InputStream fileStream = new FileInputStream(file);
			if(file.endsWith(".gz"))
				s = new BufferedReader(new InputStreamReader(new GZIPInputStream(fileStream)));
			else
				s = new BufferedReader(new InputStreamReader(fileStream));
		}catch(IOException e)
		{
			e.printStackTrace();
		}
	}
	
	public void close()
	{
		try {
			s.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	public boolean hasMore()
	{
		try {
			return s.ready();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return false;
	}
	
	public String getLine()
	{
		try {
			return s.readLine();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return null;
	}
}
