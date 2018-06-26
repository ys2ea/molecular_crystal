import java.util.*;
import java.io.*;
import static java.lang.Math.*;

//A code to seperate atoms into molecules.
//Takes account of PBC.
//Bond length criterial is set to be 1.5A in all cases, better checking should be added if needed


class GroupMol {
	List<Molecule> mollist;
	int nmol;
	int natom;
	double lx;
	double ly;
	double lz;
	static double blength = 1.5;
	
	public GroupMol(String coordinname, double lx, double ly, double lz) {
		String line = "";
		nmol = 0;
		natom = 0;
        	mollist = new ArrayList<Molecule>();
        	//here the mollist is not a molecule in the middle of the run, but is a list of fragments that keep merging the elements to get final state
        	
		try {
			BufferedReader br = new BufferedReader(new FileReader(coordinname));
		
			while ((line = br.readLine()) != null) {
				natom ++;
				Scanner scan = new Scanner(line);
				
				boolean create_new_frag = false;
				
				String t = scan.next();
				double x = scan.nextDouble();
				double y = scan.nextDouble();
				double z = scan.nextDouble();
				double dx, dy, dz, dist;
			
				//System.out.println("( " + t + ", " + natom + ", " + x + ", " + y + ", " + z + " )");
				Atom a = new Atom(t, x, y, z);
				
				List<Integer> bondlist = new ArrayList<Integer>();		//list of fragments an atom is bonded with, if more than one, merge them
				 											//it's important that the frag index is in increasing order, to ensure correct removal
				if(nmol == 0) {
					create_new_frag = true;
				}
				
				else {
					for(int i = 0; i < mollist.size(); i++) {
						for(Atom b : (mollist.get(i).atomlist) ) {
							dx = x - b.x;
							dy = y - b.y;
							dz = z - b.z;
							
							if(dx > lx/2.) dx -= lx;
							else if(dx < -lx/2.) dx += lx;
							if(dy > ly/2.) dy -= ly;
							else if(dy < -ly/2.) dy += ly;
							if(dz > lz/2.) dz -= lz;
							else if(dz < -lz/2.) dz += lz;
							
							dist = sqrt(dx*dx + dy*dy + dz*dz);
							if(dist < blength) {
								bondlist.add(i);
								break;
							}
						}
					}
					
					if(bondlist.size() == 0) create_new_frag = true;
				}
				
				if(create_new_frag) {
					nmol ++;
					Molecule m = new Molecule(lx,ly,lz);
					m.addAtom(a);
					//System.out.println("NEW!!");
					mollist.add(m);
				}
				
				//merge fragments and add current atom to the end
				else {
					if(bondlist.size()>1) {
					//System.out.println("bond list: " + bondlist);
						for(int i = 1; i < bondlist.size(); i ++) {
							mollist.get(bondlist.get(0)).atomlist.addAll(mollist.get(bondlist.get(i)).atomlist);
						}
					
						for(int j = bondlist.size()-1; j > 0; j --) {
							int rmid = bondlist.get(j);     // ----------------------------------------error when not defining rmid, why/?????????????????????????????????
							
							//System.out.println("mol # was: " + mollist.size() + " removing " + bondlist.get(j));
							mollist.remove(rmid);
							//System.out.println("mol # reduces to: " + mollist.size());
						}
						
					}
					mollist.get(bondlist.get(0)).addAtom(a);
				}
				
				//int atotal = 0;
				//for(Molecule m : mollist) atotal += m.natom;
				//System.out.println("Current atom #: " + atotal + " Mol #: " + mollist.size());
			}
			br.close();
			
			for(Molecule m : mollist) {m.print_info();}
			
		} catch  (FileNotFoundException e) {      
			System.err.println("Unable to find the file: fileName");
        	} catch (IOException e) {
            	System.err.println("Unable to read the file: fileName");
        	}
        	
        	
        	
	}
	

	public static void main(String[] argc) {
		GroupMol naph = new GroupMol("solid.xyz", 7.474, 6.2089, 37.2999);
		//GroupMol naph1 = new GroupMol("water.xyz", 30, 30, 37.2999);
	}
}
