#include "XYbilayer.h"
#include <string.h>
#include <alps/alea.h>
#include <iostream>
#include <fstream>
#include <math.h>

using namespace boost;


using namespace std;
  
bilayer::bilayer(const alps::ProcessList& where,const alps::Parameters& p,int node) : alps::scheduler::LatticeMCRun<graph_type>(where,p,node),
    L(static_cast<alps::uint32_t>(p["L"])),      // Linear lattice size
    IncStep(static_cast<string>(p["IncStep"])),
    IncNo(static_cast<alps::uint32_t>(p["IncNo"])),
    N(2*L*L),
    Nb(5*L*L),
    therm(static_cast<alps::uint64_t>(p["THERMALIZATION"])),
    beta(static_cast<double>(p["beta"])),
    J(static_cast<double>(p["Interlayer"])),
    af(true),
    sweeps(static_cast<alps::uint64_t>(p["SWEEPS"])),    // # of simulation steps
    sweeps_done(0),
    M1(static_cast<alps::uint32_t>(p["M"])),
    M2(static_cast<alps::uint32_t>(p["M"])),
    n1(0),
    n2(0),
    EnsGlued(true)
 {
    cout<<"start init"<<" "<<L<<endl;
    geom.resize(N,false);
    edge.resize(N,false);
    sit=sites().first;
    for (int i=0; i<IncStep.length(); ++i,++sit) {
        if (IncStep[i]=='1') {
            geom[*sit]=true;
            ++sit;
            geom[*sit]=true;
        }
        else if (IncStep[i]=='2') {
            edge[*sit]=true;
            ++sit;
            edge[*sit]=true;
        }
        else {
            ++sit;
        }
    }

    opstring.resize(M1+M2,0);
    bondstring.resize(M1+M2);
    BotCon.resize(N);
    maxel=4*(M1+M2+N);
   	vertexlist.resize(maxel);
    spin.resize(N,false);
    spin2.resize(N,false);
	for (int i=0; i<N; ++i) {
		spin[i]=(bool) random_int(2);
		spin2[i]=( (IsInA(i)) ? spin[i] : random_int(2) );
    }
    measurements << alps::RealObservable("ED"); //Ensemble divided at edge was S1
    measurements << alps::RealObservable("EG"); //Ensemble glued at edge was S2
}

void bilayer::fillgeo(int e1, int orient, std::vector<bool> *arr) {
	//int e1= mapx[i]+L*mapy[i];
	(*arr)[e1]=true;
	(*arr)[e1+1]=true;
	switch (orient) {
		case 0:
			(*arr)[e1+2]=true;
			(*arr)[e1+3]=true;
            /** 2x2x2 elements **/
			(*arr)[e1+2*L]=true;
			(*arr)[e1+1+2*L]=true;
			(*arr)[e1+2+2*L]=true;
			(*arr)[e1+3+2*L]=true;
            /**/
			break;
		case 1:
			(*arr)[e1+2]=true;
			(*arr)[e1+3]=true;
			(*arr)[e1+4]=true;
			(*arr)[e1+5]=true;
			break;
		case 2:
			(*arr)[e1+2*L]=true;
			(*arr)[e1+2*L+1]=true;
			(*arr)[e1+4*L]=true;
			(*arr)[e1+4*L+1]=true;
	}
}

bool bilayer::is_thermalized() const {
        return (sweeps_done>=therm);
}

double bilayer::work_done() const {
        return (is_thermalized() ? (sweeps_done-therm)/double(sweeps) :0.);
}


void bilayer::diagonal_update() {
	bond_descriptor b;
    //BotCon.clear();
    //vertexlist.clear();
	for (int i=0; i<N; ++i) {
		BotCon[i]=2*i;
	}

    for (int p=0; p<M1; ++p) {
        //cout<<opstring[p]<<endl;
        if (opstring[p]==2) {
            b=bondstring[p];
            spin[source(b)]=!spin[source(b)];
            spin[target(b)]=!spin[target(b)];
            base=4*p+4*N;
            vertexlist[base]=BotCon[source(b)];
            vertexlist[base+1]=BotCon[target(b)];
            vertexlist[BotCon[source(b)]]=base;
            vertexlist[BotCon[target(b)]]=base+1;
			BotCon[source(b)]=base+2;
			BotCon[target(b)]=base+3;
            continue;
        }
		if (opstring[p]%2==1) { //1 or 3
		    b=bondstring[p];
            fact = (bond_type(b))? J : 1.0 ;            //bond_type = 1 for interlayer bond, = 0 else
			if (random_01()<(M1-n1+1)/(Nb*beta*0.5*fact )) { // 0.25 for XY interaction
				--n1;
				opstring[p]=0;
			}
            else {
                base=4*p+4*N;
                vertexlist[base]=BotCon[source(b)];
                vertexlist[base+1]=BotCon[target(b)];
                vertexlist[BotCon[source(b)]]=base;
                vertexlist[BotCon[target(b)]]=base+1;
                BotCon[source(b)]=base+2;
                BotCon[target(b)]=base+3;
            }
      
            continue;
		}
		b=bond(random_int(Nb));
        /*
		if (spin[source(b)]==spin[target(b)])
			continue; */ //XY
        fact = (bond_type(b))? Nb*beta*0.5*J/(M1-n1) : Nb*beta*0.5/(M1-n1) ;            //bond_type = 1 for interlayer bond, = 0 else
        //cout<<bond_type(b)<<endl;
		if ((opstring[p]==0)&&((fact>1)||(random_01()<fact))) {
			opstring[p]=(spin[source(b)]==spin[target(b)])? 3 : 1;
            bondstring[p]=b;
            ++n1;
            base=4*p+4*N;
            vertexlist[base]=BotCon[source(b)];
            vertexlist[base+1]=BotCon[target(b)];
            vertexlist[BotCon[source(b)]]=base;
            vertexlist[BotCon[target(b)]]=base+1;
            BotCon[source(b)]=base+2;
            BotCon[target(b)]=base+3;
        }
    }

    for (int i=0; i<N; ++i) {
        if (!IsInA(i)) {
            if ((spin[i]!=spin2[i])) {
                spin[i]=!spin[i];
                spin2[i]=!spin2[i];
            }
            vertexlist[BotCon[i]]=2*i+1;
            vertexlist[2*i+1]=BotCon[i];
        }
        else {
            vertexlist[BotCon[i]]=2*N+2*i+1;
            vertexlist[2*N+2*i+1]=BotCon[i];
        }
        BotCon[i]=2*N+2*i;
    }

    for (int p=M1; p<M1+M2; ++p) {
        if (opstring[p]==2) {
		    b=bondstring[p];
            spin[source(b)]=!spin[source(b)];
            spin[target(b)]=!spin[target(b)];
            base=4*p+4*N;
            vertexlist[base]=BotCon[source(b)];
            vertexlist[base+1]=BotCon[target(b)];
            vertexlist[BotCon[source(b)]]=base;
            vertexlist[BotCon[target(b)]]=base+1;
			BotCon[source(b)]=base+2;
			BotCon[target(b)]=base+3;
            continue;
        }
		if (opstring[p]%2==1) {
		    b=bondstring[p];
            fact = (bond_type(b))? J : 1.0 ;
			if (random_01()<(M2-n2+1)/(Nb*beta*0.5*fact )) { //XY bilayer
				--n2;
				opstring[p]=0;
			}
            else {
                base=4*p+4*N;
                vertexlist[base]=BotCon[source(b)];
                vertexlist[base+1]=BotCon[target(b)];
                vertexlist[BotCon[source(b)]]=base;
                vertexlist[BotCon[target(b)]]=base+1;
                BotCon[source(b)]=base+2;
                BotCon[target(b)]=base+3;
            }
     
            continue;
		}
		b=bond(random_int(Nb));
        /*
		if (spin[source(b)]==spin[target(b)])
			continue; */ //XY
        fact = (bond_type(b))? Nb*beta*0.5*J/(M2-n2) : Nb*beta*0.5/(M2-n2) ; 
		if ((opstring[p]==0)&&((fact>1)||(random_01()<fact)))  {
			opstring[p]=(spin[source(b)]==spin[target(b)])? 3 : 1;
            bondstring[p]=b;
            ++n2;
            base=4*p+4*N;
            vertexlist[base]=BotCon[source(b)];
            vertexlist[base+1]=BotCon[target(b)];
            vertexlist[BotCon[source(b)]]=base;
            vertexlist[BotCon[target(b)]]=base+1;
            BotCon[source(b)]=base+2;
            BotCon[target(b)]=base+3;
        }
    }

    for (int i=0; i<N; ++i) {
        if (!IsInA(i)) {
            vertexlist[BotCon[i]]=2*N+2*i+1;
            vertexlist[2*N+2*i+1]=BotCon[i];
            if ((spin[i]!=spin2[i])) {
                spin[i]=!spin[i];
                spin2[i]=!spin2[i];
            }
        }
        else {
            vertexlist[BotCon[i]]=2*i+1;
            vertexlist[2*i+1]=BotCon[i];
        }
    }	
}


void bilayer::dostep() {
	sweeps_done++;
    diagonal_update();

	seen.clear();
    seen.resize(maxel,0);
	fliploops();
	if (isboth()) {
		EnsGlued=!EnsGlued;
        if (!EnsGlued) {
            for (int i=0; i<N; ++i) {
                if (edge[i])
                    spin2[i]=spin[i];
            }
        }
	}
    if (is_thermalized()) {
        do_measurement();
    }
    else {
        if ((M1-n1)<n1/10) {
            opstring.insert(opstring.begin()+M1,n1/10,0);
            bondstring.insert(bondstring.begin()+M1,n1/10,bond(0));
            M1+=n1/10;
            maxel=4*(M1+M2+N);
            vertexlist.resize(maxel);
            cerr << J << ": M1 too small, now "<< M1 << ", sweep=" << sweeps_done << endl;

        }
        if ((M2-n2)<n2/10) {
            M2+=n2/10;
            opstring.resize(M1+M2,0);
            bondstring.resize(M1+M2);
            maxel=4*(M1+M2+N);
            vertexlist.resize(maxel);
            cerr << J << ": M2 too small, now "<< M2 << endl;
        } 
    }
}

void bilayer::check_list() {
    for (int i=-4*(M2+N); i<4*(M1+N); i++) {
        if ((vertexlist[i])&&(vertexlist[vertexlist[i]]!=i))
            cout<<"Problem bei "<<i<<endl;
    }
    for (int i=-2*N; i<2*N; i++)
        if (!vertexlist[i])
            cout<<"Problem, fehlt "<<i<<endl;
}

void bilayer::do_measurement() {
    measurements["ED"] << double(!EnsGlued);
    measurements["EG"] << double(EnsGlued);
}

bool bilayer::IsInA(int x) {
	return ((edge[x])&&(EnsGlued))||(geom[x]);
}


void bilayer::go_loop(int p, bool flip) {
	seen[p]=1;
	if (p<2*N) {
        sn=(p%(2*N))/2;
        spin[sn]=flip!=spin[sn];
		seen[p^1]=1;
		if (!seen[vertexlist[p^1]])
		go_loop(vertexlist[p^1],flip);
	}
    else if (p<4*N) {
        sn=(p%(2*N))/2;
       	spin2[sn]=flip!=spin2[sn];
		seen[p^1]=1;
		if (!seen[vertexlist[p^1]])
		go_loop(vertexlist[p^1],flip);
    }
	else {
        switch (opstring[(p-4*N)/4]) {
            case 1: 
                if (flip) 
                    opstring[(p-4*N)/4]=2; // 1 <--> 2
                seen[p^1]=1;
                if (!seen[vertexlist[p^1]])
                    go_loop(vertexlist[p^1],flip);
                /*
                if (!seen[p^2])
                    go_loop(p^2, random_int(2));*/
                break;
            case 3:
                if (flip) 
                    opstring[(p-4*N)/4]=2;
                seen[p^3]=1;
                if (!seen[vertexlist[p^3]])
                    go_loop(vertexlist[p^3],flip);
                /*
                if (!seen[p^2])
                    go_loop(p^2, random_int(2));*/
                break;
            case 2:
                take_af=random_int(2);
                if (seen[p^1])
                    take_af=false;
                if (seen[p^3])
                    take_af=true;
                if (take_af) { //decide with prob=0.5 to which diag op to return
                    if (flip) 
                        opstring[(p-4*N)/4]=1;
                    seen[p^1]=1;
                    if (!seen[vertexlist[p^1]])
                        go_loop(vertexlist[p^1],flip);
                    break;
                }
                else {
                    if (flip) 
                        opstring[(p-4*N)/4]=3;
                    seen[p^3]=1;
                    if (!seen[vertexlist[p^3]])
                        go_loop(vertexlist[p^3],flip);
                    break;
                }
        }
    }
}


void bilayer::fliploops() {
    /*
	for (int i=0; i<4*N; i+=2)
		if (!seen[i])
			go_loop(i,random_int(2));
    */
    go_loop(random_int(4*(M1+M2+N)),true);
}

bool bilayer::isboth() {
	for (int i=0; i<N; ++i) 
		if (edge[i]&&(spin2[i]!=spin[i]))
			return false;
	return true;}


void bilayer::save(alps::ODump& dump) const
{
    std::vector<int> intbondstring(M1+M2);
    for (int i=0; i<M1+M2; ++i) 
        intbondstring[i]=edge_index(bondstring[i]);
    dump<<sweeps_done<<M1<<M2<<n1<<n2<<opstring<<intbondstring<<EnsGlued<<spin<<spin2;
}

void bilayer::load(alps::IDump& dump) 
{
    std::vector<int> intbondstring;
    dump>>sweeps_done>>M1>>M2>>n1>>n2>>opstring>>intbondstring>>EnsGlued>>spin>>spin2;
    bondstring.resize(M1+M2);
    for (int i=0; i<M1+M2; ++i) 
        bondstring[i]=bond(intbondstring[i]);
    BotCon.resize(N);
   	vertexlist.resize(4*(M1+M2+N));
    maxel=4*(M1+M2+N);
    //cout<<M1<<" "<<sweeps_done<<" "<<N<<endl;

}

void bilayer::print_copyright(std::ostream & out)
{
        out << " copyright (c) by Johannes Helmes\n";
}
