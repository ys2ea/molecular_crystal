import java.util.*;
import java.io.*;
import static java.lang.Math.*;

//A code to seperate atoms into molecules.
//Takes account of PBC.
//Bond length criterial is set to be 1.5A in all cases, better checking should be added if needed
//Now it also unwraps a molecule. It prepars molecule and charges.

class GroupMol {
	List<Molecule> mollist;
	int nmol;
	int natom;
	double lx;
	double ly;
	double lz;
	static double blength = 1.5;
	//symmetry defines how to find coord of all molecules in a unit cell, from a give ref molecule
	//the coord of ith mol is given by (x*fxi + dxi, y*fyi + dyi, z*fzi + dzi)
	//sym1[i][d] = fdi, sym2[i][d] = ddi
	static double[][] symmetry1 = {{1,1,1},{-1,1,1},{1,-1,1},{1,1,-1}};
	static double[][] symmetry2 = {{0,0,0},{0.5,-0.5,0},{0,1.5,-0.5},{0.5,0,0.5}};
	public GroupMol(String coordinname, double llx, double lly, double llz) {
		String line = "";
		nmol = 0;
		natom = 0;
		lx = llx;
		ly = lly;
		lz = llz;
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
				double c = 0;
				if(scan.hasNextDouble()) c = scan.nextDouble();
				double dx, dy, dz, dist;
			
				//System.out.println("( " + t + ", " + natom + ", " + x + ", " + y + ", " + z + " )");
				Atom a = new Atom(t, natom, x, y, z, c);
				
				List<Integer> bondlist = new ArrayList<Integer>();		//list of fragments an atom is bonded with, if more than one, merge them
				 											//it's important that the frag index is in increasing order, to ensure correct removal
				if(nmol == 0) {
					create_new_frag = true;
				}
				
				else {
					for(int i = 0; i < mollist.size(); i++) {
						boolean connect = false;
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
								if(!connect) {
									bondlist.add(i);
									connect = true; }
								b.add_nb(a);
								a.add_nb(b);
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
			
			
		} catch  (FileNotFoundException e) {      
			System.err.println("Unable to find the file: fileName");
        	} catch (IOException e) {
            	System.err.println("Unable to read the file: fileName");
        	}
        	
        	
        	
	}
	
	//qmlist is list of molecules that are treated with qm.
	//There are 4 molecules in a unit cell. A molecule is noted by (nx, ny, nz, id), id= 0, 1, 2, 3. 
	public void print_cfile(String name_qm, String name_c, List<int[]> qmlist, int xmin, int xmax, int ymin, int ymax, int zmin, int zmax) throws IOException {
		FileWriter fw1=new FileWriter(name_qm);
		FileWriter fw2=new FileWriter(name_c);
		
		//for(int t = 0; t < qmlist.size(); t ++) if(qmlist.get(t).length!=4) System.out.println("QMLIST SIZE WRONG!!!!!!!!");
		
		String line;
		Molecule m = mollist.get(0);
		
		for(int k = xmin; k <= xmax ; k ++) {
			for(int n = ymin; n <= ymax ; n ++) {
				for(int l = zmin; l <= zmax ; l ++) {
					for(int mid = 0; mid < 4; mid ++) {
						boolean isqm = false;
						//decide whether a molecule is qm part.
						for(int idq = 0; idq < qmlist.size(); idq ++) {
							if(k==qmlist.get(idq)[0] && n==qmlist.get(idq)[1] && l==qmlist.get(idq)[2] && mid==qmlist.get(idq)[3]) {
								isqm = true;
								System.out.println("Q id: " + k + ", " + n + ", " + l + ", " + mid);
								break;
							}
						}
					
						if(!isqm) {
							for(int i =0; i < m.atomlist.size(); i++) {
								line = String.format("%s  %f  %f  %f  %f  %n", m.atomlist.get(i).type, (symmetry1[mid][0]*m.atomlist.get(i).x+symmetry2[mid][0]*lx+lx*k), 
								 (symmetry1[mid][1]*m.atomlist.get(i).y+symmetry2[mid][1]*ly+ly*n), (symmetry1[mid][2]*m.atomlist.get(i).z+symmetry2[mid][2]*lz+lz*l), m.atomlist.get(i).charge);
								fw2.write(line);
							}
						}
					
						else {
							for(int i =0; i < m.atomlist.size(); i++) {
								line = String.format("%s  %f  %f  %f  %n", m.atomlist.get(i).type, (symmetry1[mid][0]*m.atomlist.get(i).x+symmetry2[mid][0]*lx+lx*k), 
								 (symmetry1[mid][1]*m.atomlist.get(i).y+symmetry2[mid][1]*ly+ly*n), (symmetry1[mid][2]*m.atomlist.get(i).z+symmetry2[mid][2]*lz+lz*l));
								fw1.write(line);
							}
						}
					} //end mid
				}  //end l
			} //end n
		}// end k

		fw1.close();
		fw2.close();
	}
	
	public static void main(String[] argc) throws IOException {
		GroupMol naph = new GroupMol("mol.xyz", 7.474, 6.2089, 37.2999);
		//System.out.println("Total mol: " + naph.mollist.size());
		//for(Molecule m : naph.mollist) {m.print_info();}
		naph.mollist.get(0).unwrap();
		//naph.mollist.get(0).print_info();
		int[] qm1 = {0,0,0,1};
		int[] qm2 = {0,0,0,0};
		List<int[]> qmlist = new ArrayList<int[]>();
		qmlist.add(qm1);
		qmlist.add(qm2);
		naph.print_cfile("qm_coord.xyz", "charge_coord.xyz", qmlist, -1, 1, -1, 1, -1, 1);
	}
}
