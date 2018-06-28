import java.util.*;
class Atom {
	String type;
	int idx;
	double charge;
	double x;
	double y;
	double z;
	boolean tag;		//an extra tag to play with;
	List<Atom> nb;
	
	public Atom(String t, int i, double xin, double yin, double zin, double c) {
		type = t;
		idx = i;
		x = xin;
		y = yin;
		z = zin;
		charge = c;
		nb = new ArrayList<Atom>();
		tag = false;
	}
	
	public void set_coord(double xnew, double ynew, double znew) {
		x = xnew;
		y = ynew;
		z = znew;
	}
	
	public void add_nb(Atom b) {
		nb.add(b);
	}
}
