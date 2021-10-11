#ifndef _POINT_H_
#define _POINT_H_
#include<vector>
#include<algorithm>
#include<CGAL/basic.h>

struct Point {
	typedef double Cor;
	typedef double Val;
	typedef double Rw;
	typedef std::vector<Cor>::reference Ref;
	typedef std::vector<Cor>::const_reference Cref;
	typedef std::vector<Cor>::iterator It;
	typedef std::vector<Cor>::const_iterator Cit;
	
	std::vector<Cor> dcomp;
	Val dval=0;
	Rw drw=1;

	Point():
		dcomp(),
		dval(0),
		drw(1)
	{};
	Point(const Point& p):
		dcomp(p.begin(),p.end()),
		dval(p.val()),
		drw(p.rw())
	{};
	template <class inputIt> Point(inputIt beg,inputIt past):
		dcomp(beg,past),
		dval(0),
		drw(1)
	{}
	template <class inputIt> Point(inputIt beg, inputIt past, double & val):
		dcomp(beg,past),
		dval(val),
		drw(1)
	{}
	Point(std::size_t const & n, Cor c):
		dcomp(n,c),
		dval(0),
		drw(1)
	{};

	Ref operator[] (std::size_t n) {return dcomp[n];}
	Cref operator[] (std::size_t n) const {return dcomp[n];}
	It begin() {return dcomp.begin();};
	Cit begin() const {return dcomp.begin();};
	It end() {return dcomp.end();};
	Cit end() const {return dcomp.end();};
	void push_back(double & val) {dcomp.push_back(val);};
	void val(Val const & in) {dval=in;};
	Val const & val() const {return dval;};
	void rw(Rw const & in) {drw=in;};
	Rw const & rw() const {return drw;};
	std::size_t const dims() const {return dcomp.size();};

	bool operator==(Point const& p) const {
		return (std::equal(begin(),end(),p.begin()) && val()==p.val() && rw()==p.rw());
	}
	bool operator!=(Point const& p) const {
		return !(*this == p);
	}
};

namespace CGAL {
  template <> struct Kernel_traits<Point> {
    struct Kernel {
      typedef double FT;
      typedef double RT;
    };
  };
}


struct Construct_coord_iterator {
  typedef  Point::Cit result_type;
  Point::Cit operator()(const Point& p) const
  { return p.begin(); }

  Point::Cit operator()(const Point& p, int)  const
  { return p.end(); }
};
#endif
