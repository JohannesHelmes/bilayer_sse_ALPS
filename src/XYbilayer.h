#ifndef HBERG
#define HBERG
#include <boost/random.hpp>
#include <alps/scheduler.h>
#include <alps/lattice/graph_helper.h>
#include <alps/alea.h>
#include <alps/scheduler/montecarlo.h>
#include <map>
#include <string>

#include <boost/unordered_set.hpp>

typedef alps::scheduler::LatticeMCRun<>::graph_type graph_type;
    
using namespace std;
        
class bilayer : public alps::scheduler::LatticeMCRun<graph_type>{

public:
    bilayer(const alps::ProcessList& where,const alps::Parameters& p,int node);
    static void print_copyright(std::ostream &);
    void save(alps::ODump& dump) const;
    void load(alps::IDump& dump);
    void dostep();
    bool is_thermalized() const;
    double work_done() const;

private:
    void do_measurement();
    void check_list();
    void diagonal_update();


	std::vector<int> opstring;
    std::vector<bond_descriptor> bondstring;
	std::vector<int> BotCon;
	std::vector<bool> spin;
	std::vector<bool> spin2;
	std::vector<bool> geom;
	std::vector<bool> edge;
	alps::uint64_t sweeps,sweeps_done,therm;
	int IncNo,L,La,N,n1,n2,M1,M2,Nb,blength,maxel;
	int edge1,base,sn;
	bool af,EnsGlued,even,take_af;
	std::vector<int> seen;
	std::vector<int> vertexlist;
	double beta,J,fact;
    string IncStep;
    site_iterator sit;

	void fillgeo(int e1, int orient, std::vector<bool> *arr);
	void fliploops();
	void go_loop(int p, bool flip);
	bool IsInA(int x);
	bool isboth();

};
typedef alps::scheduler::SimpleMCFactory<bilayer> bilayerFactory;
#endif /* HBERG */
