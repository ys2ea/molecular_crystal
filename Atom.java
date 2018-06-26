class Atom {
	String type;
	
	double x;
	double y;
	double z;
	
	public Atom(String t, double xin, double yin, double zin) {
		type = t;
		x = xin;
		y = yin;
		z = zin;
	}
	
	public void set_coord(double xnew, double ynew, double znew) {
		x = xnew;
		y = ynew;
		z = znew;
	}
}
