import java.util.*;
import java.io.*;
import static java.lang.Math.*;

class Molecule {
	List<Atom> atomlist;
	int natom;
	double lx;
	double ly;
	double lz;
	static double bl = 1.5;
	public Molecule(double x, double y, double z) {
		lx = x;
		ly = y;
		lz = z;
		natom = 0;
		atomlist = new ArrayList<Atom>();
	}
	
	public void addAtom(Atom a) {
		natom ++;
		atomlist.add(a);
	}
	
	public void print_info() {
		System.out.printf("# of Atom: %d  list: ", atomlist.size());
		double totalx = 0., totaly = 0., totalz = 0.; 
		for(int i =0; i < atomlist.size(); i++) {
			System.out.printf("(%s,%d,%d)", atomlist.get(i).type, atomlist.get(i).idx, atomlist.get(i).nb.size());
			totalx += atomlist.get(i).x; totaly += atomlist.get(i).y; totalz += atomlist.get(i).z;
		}
		System.out.printf("\n");
		double x1=0., y1=0., z1=0., x2=0., y2=0., z2=0.;
		boolean first = false;
		System.out.println("total coord: " + totalx/lx/atomlist.size() + ", " + totaly/ly/atomlist.size() + ", " + totalz/lz/atomlist.size());
		for(int i =0; i < atomlist.size(); i++) {
		
			if(atomlist.get(i).type.charAt(0) == 'O') 
			{
				System.out.println(atomlist.get(i).type + " coord: (" + atomlist.get(i).x + ", " + atomlist.get(i).y + ", " + atomlist.get(i).z + ")");
				if(!first) {
					x1 = atomlist.get(i).x; y1 = atomlist.get(i).y; z1 = atomlist.get(i).z; 
					first = true;
				}
				
				else { x2 = atomlist.get(i).x; y2 = atomlist.get(i).y; z2 = atomlist.get(i).z; }
			}
		}
		
		
		//dxo is from o2 to o1
		
		double dxo = -0.008114337563780483;
		double dyo = -0.3158010254036201;
		double dzo = -0.9487907408274614;
		
	//	double dxo = x1-x2;
	//	double dyo = y1-y2;
	//	double dzo = z1-z2;
		
	//	if(dxo > lx/2.) dxo -= lx;
	//	else if(dxo < -lx/2.) dxo += lx;
	//	if(dyo > ly/2.) dyo -= ly;
	//	else if(dyo < -ly/2.) dyo += ly;
	//	if(dzo > lz/2.) dzo -= lz;
	//	else if(dzo < -lz/2.) dzo += lz;

		
	//	double mdo = sqrt(dxo*dxo+dyo*dyo+dzo*dzo);
	//	dxo /= mdo;
	//	dyo /= mdo;
	//	dzo /= mdo;
		
	//	System.out.println("Z vector: (" + dxo + ", " + dyo + ", " + dzo + ")");
		
		//extra work to find atoms bonded to O
		for(int i =0; i < atomlist.size(); i++) {
			double dx = atomlist.get(i).x - x1;
			double dy = atomlist.get(i).y - y1;
			double dz = atomlist.get(i).z - z1;
			
			if(dx > lx/2.) dx -= lx;
			else if(dx < -lx/2.) dx += lx;
			if(dy > ly/2.) dy -= ly;
			else if(dy < -ly/2.) dy += ly;
			if(dz > lz/2.) dz -= lz;
			else if(dz < -lz/2.) dz += lz;

			if(sqrt(dx*dx+dy*dy+dz*dz)<bl && atomlist.get(i).type.charAt(0) != 'O') {
				System.out.println("Bonded to O1: " + atomlist.get(i).type + " dx: " + dx + " dy: " + dy + " dz: " + dz);
				double r = dx*dxo + dy*dyo + dz*dzo;
				dx -= dxo*r;
				dy -= dyo*r;
				dz -= dzo*r;
				double nd = sqrt(dx*dx + dy*dy+dz*dz);
				System.out.println("projected vector: dx: " + dx/nd + " dy: " + dy/nd + " dz: " + dz/nd + " product: " + (dx*dxo+dy*dyo+dz*dzo) );
			}
			
			else {
				dx = atomlist.get(i).x - x2;
				dy = atomlist.get(i).y - y2;
				dz = atomlist.get(i).z - z2;
				
				if(dx > lx/2.) dx -= lx;
				else if(dx < -lx/2.) dx += lx;
				if(dy > ly/2.) dy -= ly;
				else if(dy < -ly/2.) dy += ly;
				if(dz > lz/2.) dz -= lz;
				else if(dz < -lz/2.) dz += lz;
				
				if(sqrt(lx*lx+ly*ly+lz*lz)<bl && atomlist.get(i).type.charAt(0) != 'O') {
					System.out.println("Bonded to O2: " + atomlist.get(i).type + " dx: " + dx + " dy: " + dy + " dz: " + dz);
					double r = dx*dxo + dy*dyo + dz*dzo;
					dx += dxo*r;
					dy += dyo*r;
					dz += dzo*r;
					double nd = sqrt(dx*dx + dy*dy+dz*dz);
					System.out.println("projected vector: dx: " + dx/nd + " dy: " + dy/nd + " dz: " + dz/nd + " product: " + (dx*dxo+dy*dyo+dz*dzo));					
				}
			}
		}
			
	}
	
	public void print_to_file(String name) throws IOException {
		FileWriter fw=new FileWriter(name);
		String line;
		for(int i =0; i < atomlist.size(); i++) {
			line = String.format("%s   %f  %f  %f  %f %n", atomlist.get(i).type, atomlist.get(i).x, atomlist.get(i).y, atomlist.get(i).z, atomlist.get(i).charge);
			fw.write(line);
		}
		fw.close();
	}
	
	public void reset_tag() {
		for(Atom a : atomlist) a.tag = false;
	}
	
	//unwrap a molecule from PBC	
	public void unwrap() {
		if(lx<bl*3. || ly<bl*3. || lz<bl*3.) { System.out.println("System too small to unwrap!!"); System.exit(1); }
		
		reset_tag();
		List<Atom> stack = new ArrayList<Atom>();
		double dx, dy, dz;
		
		stack.add(atomlist.get(0));
		atomlist.get(0).tag = true;
		while(stack.size()!=0 && stack.size() < 100) {
			int last = stack.size() - 1;

			boolean rm = true;
			for(Atom b : stack.get(last).nb) {
				if(!b.tag) {
					stack.add(b); 
					b.tag = true;
					rm = false;
					dx = stack.get(last).x - b.x;
					dy = stack.get(last).y - b.y;
					dz = stack.get(last).z - b.z;
				
					if(dx > lx-bl) b.x += lx;
					else if(dx < bl-lx) b.x -= lx;
					if(dy > ly-bl) b.y += ly;
					else if(dy < bl-ly) b.y -= ly;
					if(dz > lz-bl) b.z += lz;
					else if(dz < bl-lz) b.z -= lz;
				}
			}
			if(rm) stack.remove(last);
		}
	}
}
