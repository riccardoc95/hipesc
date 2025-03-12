#include <iostream>
#include <unordered_map>
#include "robin_hood.h"
#include <algorithm>
#include <mutex>
#include <limits>
#include <stdint.h>
#include <unordered_map>
#include <string>
#include <vector>
#include <set>
#include <utility>
#include <fstream>
#include <stdio.h>
#include <atomic>
#include <math.h>
#include <cmath>
#include <cstdint>

#include "utils.h"

#include "spoa/spoa.hpp"
#include "ssw_cpp.h"

typedef uint32_t kmer;

struct localisation {
    uint32_t read_id;
    int32_t position;
};

//TODO REMOVE HASHTABLES

using namespace std;
using uint = unsigned int;


using namespace std;
typedef unsigned int uint;



kmer str2num(const string& str){
    kmer res(0);
    for(uint i(0);i<str.size();i++){
        res<<=2;
        switch (str[i]){
            case 'A':res+=0;break;
            case 'C':res+=1;break;
            case 'G':res+=2;break;
            default:res+=3;break;
        }
    }
    return res;
}



string kmer2str(kmer k,uint32_t nuc){
    string result;
    for(uint32_t i(0);i<nuc;++i){
        switch(k%4){
            case 0:  result.push_back('A');break;
            case 1:  result.push_back('C');break;
            case 2:  result.push_back('G');break;
            case 3:  result.push_back('T');break;
        }
        k>>=2;
    }
    reverse(result.begin(),result.end());
    return result;
}


kmer nuc2int(char c){
    switch(c){
        //~ case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
    }
    return 0;
}



kmer nuc2intrc(char c){
    switch(c){
        case 'A': return 3;
        case 'C': return 2;
        case 'G': return 1;
            //~ case 'T': return 0;
    }
    return 0;
}




string intToString(uint64_t n){
    if(n<1000){
        return to_string(n);
    }
    string end(to_string(n%1000));
    if(end.size()==3){
        return intToString(n/1000)+","+end;
    }
    if(end.size()==2){
        return intToString(n/1000)+",0"+end;
    }
    return intToString(n/1000)+",00"+end;
}



void updateK(kmer& min, char nuc,kmer offsetUpdateKmer){
    min<<=2;
    min+=nuc2int(nuc);
    min%=offsetUpdateKmer;
}



char randNuc(){
    switch (rand()%4){
        case 0:
            return 'A';
        case 1:
            return 'C';
        case 2:
            return 'G';
        case 3:
            return 'T';
    }
    return 'A';
}



string rand_seq(const int n){
    string result;
    for(int32_t i(0);i<n;++i){
        result.push_back(randNuc());
    }
    return result;
}



string mutate(string& read, int n){
    for(int i(0); i<n; ++i){
        int position(rand()%read.size());
        int mutation_type(rand()%3);
        if(mutation_type==0){
            //~ if(mutation_type<4){
            read[position]=randNuc();
        }else if (mutation_type==1){
            read=read.substr(0,position)+read.substr(position+1);
        }else{
            read=read.substr(0,position)+randNuc()+read.substr(position);
        }
    }
    return read;
}


class rev_comp {
public:
    static std::string run(std::string seq);

protected:
    static rev_comp* _instance;

    static void build_instance();

    rev_comp();

    ~rev_comp();

private:
    static char complement[int('t') + 1];
};





struct score_chain {
    int32_t length;
    int32_t score;
    int32_t next_anchor;
};



typedef robin_hood::unordered_map<kmer,vector<localisation>> kmer2localisation;



void fill_index_kmers(const vector<string>& Reads,kmer2localisation& kmer_index,uint32_t kmer_size, robin_hood::unordered_map<kmer, unsigned>& merCounts, unsigned solidThresh){
    robin_hood::unordered_map<kmer, unsigned> tmpMerCounts;
    string read;
    uint32_t offsetUpdateKmer=1<<(2*kmer_size);
    robin_hood::unordered_map<kmer,bool> repeated_kmer;
    localisation here({0,0});
    for(uint32_t iR(0);iR<Reads.size();++iR){
        robin_hood::unordered_map<kmer,uint32_t> local_kmer;
        here.position=0;
        here.read_id=iR;
        read=Reads[iR];
        if(read.size()<kmer_size){continue;}
        kmer seq(str2num(read.substr(0,kmer_size)));
        kmer_index[seq].push_back(here);
        // tmpMerCounts[seq]++;
        if(++local_kmer[seq]>1){
            repeated_kmer[seq]=true;
        }
        for(uint32_t ir(0);kmer_size+ir<read.size();++ir){
            updateK(seq,read[kmer_size+ir],offsetUpdateKmer);
            ++here.position;
            kmer_index[seq].push_back(here);
            // tmpMerCounts[seq]++;
            if(++local_kmer[seq]>1){
                repeated_kmer[seq]=true;
            }
        }
    }

    for (auto p : kmer_index) {
        if (p.second.size() >= solidThresh) {
            merCounts[p.first] = p.second.size();
        }
    }

    auto it = repeated_kmer.begin();
    while(it != repeated_kmer.end()){
        kmer_index.erase(it->first);
        ++it;
    }

}



robin_hood::unordered_map<kmer,uint32_t> filter_index_kmers(kmer2localisation& kmer_index, double amount){
    robin_hood::unordered_map<kmer,uint32_t> result;
    //~ cerr<<"kmer INdex size before cleaning"<<kmer_index.size()<<endl;
    vector<kmer> to_suppress;
    vector<uint32_t> read_ids;
    auto it = kmer_index.begin();
    while(it != kmer_index.end()){
        for(uint32_t i(0);i<it->second.size();++i){
            read_ids.push_back(it->second[i].read_id);
        }
        //AVOID TO COUNT MULTIPLE OCCURENCE OF A KMER WITHIN A READ
        sort( read_ids.begin(), read_ids.end() );
        int uniqueCount =distance(read_ids.begin(), unique(read_ids.begin(), read_ids.end())) ;
        if(uniqueCount<amount){
            to_suppress.push_back(it->first);
        }else{
            result[it->first]=uniqueCount;
        }
        ++it;
        read_ids.clear();
    }
    for(uint32_t i(0);i<to_suppress.size();++i){
        kmer_index.erase(to_suppress[i]);
    }
    //~ cerr<<"kmer INdex size after cleaning"<<kmer_index.size()<<endl;
    return result;
}


void erase_current_element(vector<localisation>& V,uint n){
    for(uint i(n);i+1<V.size();++i){
        V[i]=V[i+1];
    }
    V.resize(V.size()-1);
}



void clean_suspcious_reads(kmer2localisation& kmer_index, uint read_number,double threshold){
    vector<bool> read_ok(read_number,true);
    vector<uint32_t> read_seed_number(read_number,0);
    {
        auto it = kmer_index.begin();
        while(it != kmer_index.end()){
            for(uint32_t i(0);i<it->second.size();++i){
                read_seed_number[(it->second[i].read_id)]++;
            }
            ++it;
        }
        for(uint i(0);i< read_number;++i){
            if(read_seed_number[i]<threshold){
                read_ok[i]=false;
            }
        }
    }
    //~ cout<<"goo"<<endl;
    auto it = kmer_index.begin();
    while(it != kmer_index.end()){
        for(uint32_t i(0);i<it->second.size();++i){
            if(not read_ok[it->second[i].read_id]){
                erase_current_element(it->second,i);
            }
        }
        ++it;
    }
}



bool order_according2read_id (localisation i,localisation j) { return (i.read_id<j.read_id); }



int anchors_ordered_according2reads(const kmer kmer1,const kmer kmer2,  kmer2localisation& kmer_index){
    int32_t result(0);
    auto v_loc1(kmer_index[kmer1]);
    auto v_loc2(kmer_index[kmer2]);
    //~ sort (v_loc1.begin(), v_loc1.end(), order_according2read_id);
    //~ sort (v_loc2.begin(), v_loc2.end(), order_according2read_id);
    uint32_t i1(0),i2(0);
    //BIG QUESTION HOW TO HANDLE REPEATED KMER HERE
    while(i1<v_loc1.size() and i2<v_loc2.size()){
        if(v_loc1[i1].read_id==v_loc2[i2].read_id){
            if(v_loc1[i1].position>v_loc2[i2].position){
                return -1;
                //COULD ADD A NO IF POSITIONS ARE TOO FAR LIKE IN MINIMAP
            }else{
                ++i1;
                ++i2;
                ++result;
            }
        }else if(v_loc1[i1].read_id<v_loc2[i2].read_id){
            i1++;
        }else{
            i2++;
        }
    }
    return result;
}



score_chain longest_ordered_chain_from_anchors( kmer2localisation& kmer_index, robin_hood::unordered_map<uint,score_chain>& best_chain_computed, uint32_t start, const vector<kmer>& template_read,double edge_solidity){
    if(best_chain_computed.count(start)==1){
        return best_chain_computed[start];
    }
    int32_t max_chain(-1),max_score(0);
    int32_t next_anchor(-1);
    for(uint i(start+1);i<template_read.size();++i){
        kmer next(template_read[i]);
        int score(anchors_ordered_according2reads(template_read[start],next,kmer_index));
        if(score>=edge_solidity){
            auto p=longest_ordered_chain_from_anchors(kmer_index,best_chain_computed,i,template_read,edge_solidity);
            if(p.length>max_chain){
                max_chain=p.length;
                max_score=p.score+score;
                next_anchor=i;
            }else if(p.length==max_chain and p.score+score>max_score) {
                max_score=p.score+score;
                next_anchor=i;
            }else{
            }
        }
    }

    //~ cerr<<"SCORE of "<<start<<": "<<max_chain+1<<" "<<max_score<<" "<<next_anchor<<endl;
    best_chain_computed[start]={max_chain+1,max_score,next_anchor};
    return {max_chain+1,max_score,next_anchor};
}



vector<kmer> get_template( kmer2localisation& kmer_index,const string& read,int kmer_size){
    vector<kmer> result;
    uint32_t offsetUpdateKmer=1<<(2*kmer_size);
    kmer seq(str2num(read.substr(0,kmer_size)));
    if(kmer_index.count(seq)){
        result.push_back(seq);
    }
    for(uint32_t ir(0);kmer_size+ir<read.size();++ir){
        updateK(seq,read[kmer_size+ir],offsetUpdateKmer);
        if(kmer_index.count(seq)){
            result.push_back(seq);
        }
    }
    return result;
}




vector<kmer> longest_ordered_chain( kmer2localisation& kmer_index,const vector<kmer>& template_read, double edge_solidity){
    robin_hood::unordered_map<uint,score_chain> best_chain_computed;
    vector<kmer> result;
    int32_t max_chain(0),max_score(0);
    int32_t next_anchor(-1);
    for(int32_t i(template_read.size()-1);i>=0;--i){
        auto p=longest_ordered_chain_from_anchors(kmer_index,best_chain_computed,i,template_read,edge_solidity);
        if(p.length>max_chain){
            max_chain=p.length;
            max_score=p.score;
            next_anchor=i;
        }else if(p.length==max_chain and p.score>max_score) {
            max_score=p.score;
            next_anchor=i;
        }
    }
    while(next_anchor!=-1){
        result.push_back(template_read[next_anchor]);
        next_anchor=best_chain_computed[next_anchor].next_anchor;
    }
    return result;
}



bool comparable(double x, pair<double,double> deciles){
    //~ return true;
    if(x<deciles.second+5){
        //LOW VALUE
        if(x>deciles.first-5){
            return true;
        }
        if(x/deciles.first<0.5){
            return false;
        }
        return true;
    }else{
        //High value
        if(x/deciles.second>2){
            return false;
        }
        return true;
    }
}



bool comparable(double x,double mean){
    //~ return true;
    if(abs(x-mean)<5){
        return true;
    }
    if(x/mean<0.5 or x/mean>2){
        return false;
    }
    return true;
}



//~ double mean(const vector<uint32_t>& V){
//~ double first_mean(0);
//~ uint32_t valid(0);
//~ for(uint32_t i(0);i<V.size();++i){
//~ first_mean+=V[i];
//~ }
//~ first_mean/=V.size();
//~ return first_mean;
//~ double second_mean(0);
//~ for(uint32_t i(0);i<V.size();++i){
//~ if(comparable((double)V[i],first_mean)){
//~ second_mean+=V[i];
//~ valid++;
//~ }
//~ }
//~ if(valid!=0){
//~ second_mean/=valid;
//~ }else{
//~ return first_mean;
//~ }
//~ return second_mean;
//~ }



pair<double,double> deciles( vector<uint32_t>& V){
    sort(V.begin(),V.end());
    return {V[floor((V.size()-1)*0.2)],V[ceil((V.size()-1)*0.8)]};
}



vector<double> average_distance_next_anchor(kmer2localisation& kmer_index,  vector<kmer>& anchors,robin_hood::unordered_map<kmer,uint32_t>& k_count, bool clean){
    vector<double> result;
    vector<uint32_t> v_dis;
    vector<kmer> curated_anchors;
    uint32_t min_distance(5);

    //~ {
    //~ v_dis.clear();
    //~ uint32_t sum(0),count(0);
    //~ auto v_loc1(kmer_index[anchors[i]]);
    //~ uint32_t i1(0);
    //~ while(i1<v_loc1.size()){
    //~ if(v_loc1[i1].read_id==v_loc2[i2].read_id){
    //~ v_dis.push_back(v_loc2[i1].position);
    //~ ++i1;
    //~ }
    //~ }
    //~ auto dec(deciles(v_dis));
    //~ for(uint32_t iD(0);iD<v_dis.size();++iD){
    //~ if(comparable(v_dis[iD],dec)){
    //~ sum+=v_dis[iD];
    //~ count++;
    //~ }
    //~ }
    //~ if(count==0){
    //~ cerr<<"SHOULD NOT HAPPEN"<<endl;
    //~ for(uint32_t iD(0);iD<v_dis.size();++iD){
    //~ sum+=v_dis[iD];
    //~ count++;
    //~ }
    //~ }else{
    //~ result.push_back(sum/count);
    //~ }
    //~ }

    for(uint i(0);i+1<anchors.size();++i){
        v_dis.clear();
        uint32_t sum(0),count(0);
        auto v_loc1(kmer_index[anchors[i]]);//THEY SHOULD BE READS SORTED
        auto v_loc2(kmer_index[anchors[i+1]]);
        //~ sort (v_loc1.begin(), v_loc1.end(), order_according2read_id);
        //~ sort (v_loc2.begin(), v_loc2.end(), order_according2read_id);
        uint32_t i1(0),i2(0);
        while(i1<v_loc1.size() and i2<v_loc2.size()){
            if(v_loc1[i1].read_id==v_loc2[i2].read_id){
                v_dis.push_back(v_loc2[i2].position-v_loc1[i1].position);
                ++i1;
                ++i2;
            }else if(v_loc1[i1].read_id<v_loc2[i2].read_id){
                i1++;
            }else{
                i2++;
            }
        }
        //~ if(v_dis.empty()){
        //~ }
        auto dec(deciles(v_dis));
        for(uint32_t iD(0);iD<v_dis.size();++iD){
            if(comparable(v_dis[iD],dec)){
                sum+=v_dis[iD];
                count++;
            }
        }
        if(count==0){
            // cerr<<"SHOULD NOT HAPPEN"<<endl;
            for(uint32_t iD(0);iD<v_dis.size();++iD){
                sum+=v_dis[iD];
                count++;
            }
        }else{
            result.push_back(sum/count);
        }

        //~ if(count!=0){
        //~ v_sum.push_back(sum);
        //~ v_count.push_back(count);

        //~ result.push_back(sum/count);
        //~ }else{
        //~ cerr<<"SHOULD NOT HAPPEND"<<endl;cin.get();
        //~ result.push_back(-1);
        //~ }
    }



    return result;
}



int32_t get_position(kmer2localisation& kmer_index,kmer query, uint32_t read_id){
    auto V(kmer_index[query]);
    for(uint32_t i(0);i<V.size();++i){
        if(V[i].read_id==read_id){
            return V[i].position;
        }
    }
    return -1;
}



vector<vector<string>> split_reads_old(const vector<kmer>& anchors, const vector<double>& relative_positions, const vector<string>& Reads,  kmer2localisation& kmer_index,uint32_t kmer_size){
    vector<vector<string>> result;
    for(uint32_t iR(0);iR<Reads.size();++iR){
        string read=Reads[iR];
        vector<string> split(anchors.size()+1);
        //FIRST AND LAST REGION
        int32_t anchor_position(get_position(kmer_index,anchors[0],iR));
        if(anchor_position!=-1){
            split[0]=read.substr(0,anchor_position);
        }
        anchor_position=(get_position(kmer_index,anchors[anchors.size()-1],iR));
        if(anchor_position!=-1){
            split[anchors.size()]=read.substr(anchor_position);
        }

        for(uint32_t iA(0);iA+1<anchors.size();++iA){
            int32_t anchor_position1(get_position(kmer_index,anchors[iA],iR));
            int32_t anchor_position2(get_position(kmer_index,anchors[iA+1],iR));
            if(anchor_position1!=-1){
                if(anchor_position2!=-1){
                    //REGION WITH BOtH ANCHORS
                    split[iA+1]=read.substr(anchor_position1,anchor_position2-anchor_position1);
                }else{
                    //GOT THE LEFT ANCHOR
                    //~ split[iA+1]=read.substr(anchor_position1,relative_positions[iA]);
                }
            }else{
                if(anchor_position2!=-1){
                    //GOT THE RIGHT ANCHOR
                    //~ if(anchor_position2>relative_positions[iA]){
                    //~ split[iA+1]=read.substr(anchor_position2-relative_positions[iA],relative_positions[iA]);
                    //~ }
                }
            }
        }
        result.push_back(split);
    }
    return result;
}



vector<vector<string>> split_reads(const vector<kmer>& anchors, const vector<double>& relative_positions, const vector<string>& Reads,  kmer2localisation& kmer_index,uint32_t kmer_size){
    vector<vector<string>> result(anchors.size()+1);
    if(anchors.size()==0){
        result.push_back(Reads);
        return result;
    }
    for(uint32_t iR(0);iR<Reads.size();++iR){
        //~ cerr<<endl;
        string read=Reads[iR];
        vector<string> split(anchors.size()+1);
        //FIRST AND LAST REGION
        int32_t anchor_position(get_position(kmer_index,anchors[0],iR));
        if(anchor_position!=-1){
            string chunk(read.substr(0,anchor_position));
            //~ if(abs((int)chunk.size()-relative_positions[iA])<get_position(kmer_index,anchors[0],0)*0.5){
            if(comparable(chunk.size(),get_position(kmer_index,anchors[0],0)) && chunk != ""){
                result[0].push_back(chunk);
            }
        }else{
            //~ result[0].push_back("");
        }
        anchor_position=(get_position(kmer_index,anchors[anchors.size()-1],iR));
        if(anchor_position!=-1){
            string chunk(read.substr(anchor_position));
            if(comparable(chunk.size(),Reads[0].size()-get_position(kmer_index,anchors[anchors.size()-1],0)) && chunk != ""){
                result[anchors.size()].push_back(chunk);
            }
        }else{
            //~ result[anchors.size()].push_back("");
        }
        for(uint32_t iA(0);iA+1<anchors.size();++iA){
            int32_t anchor_position1(get_position(kmer_index,anchors[iA],iR));
            int32_t anchor_position2(get_position(kmer_index,anchors[iA+1],iR));
            if(anchor_position1!=-1){
                if(anchor_position2!=-1){
                    //REGION WITH BOtH ANCHORS
                    string chunk(read.substr(anchor_position1,anchor_position2-anchor_position1));
                    //~ if(abs((int)chunk.size()-relative_positions[iA])<relative_positions[iA]*0.5){
                    if(comparable(chunk.size(), relative_positions[iA]) && chunk != ""){
                        result[iA+1].push_back(chunk);
                        //~ cerr<<chunk<<".";
                    }else{
                        //~ cerr<<"ALIEN"<<endl;
                        //~ cerr<<chunk.size()<<" "<<relative_positions[iA]<<endl;
                    }
                }else{
                    //~ cerr<<'-';
                    continue;
                    //GOT THE LEFT ANCHOR
                    string chunk(read.substr(anchor_position1,relative_positions[iA]));
                    if(comparable(chunk.size(),get_position(kmer_index,anchors[0],0)) && chunk != ""){
                        result[iA+1].push_back(chunk);
                    }else{
                        //~ cerr<<"ALIEN32"<<endl;
                    }
                }
            }else{
                if(anchor_position2!=-1){
                    //~ cerr<<'-';
                    continue;
                    //GOT THE RIGHT ANCHOR
                    if(anchor_position2>relative_positions[iA]){
                        string chunk(read.substr(anchor_position2-relative_positions[iA],relative_positions[iA]));
                        if(comparable(chunk.size(),get_position(kmer_index,anchors[0],0)) && chunk != ""){
                            result[iA+1].push_back(chunk);

                        }else{
                            //~ cerr<<"ALIEN23"<<endl;
                        }
                    }
                }else{
                    //~ cerr<<'-';
                    //~ result[iA+1].push_back("");
                }
            }
        }
    }
    return result;
}

void absoluteMAJ_consensus(vector<string>& V){
    sort(V.begin(),V.end());
    uint score(1),best_occ(0),best_score(0);
    for(uint i(0);i<V.size();++i){
        if(i+1<V.size()){
            if(V[i]!=V[i+1]){
                if(score> 0.5*V.size()){
                    V={V[i]};
                    return;
                }
                if(score>best_score){
                    best_score==score;
                    best_occ=i;
                }
                score=1;
            }else{
                score++;

            }
        }else{
            if(score> 0.5*V.size()){
                V={V[i]};
                return;
            }
        }
    }
    //~ V={V[best_occ]};
}

vector<string> consensus_SPOA( vector<string>& W, unsigned maxMSA, string path) {
    auto alignment_engine = spoa::AlignmentEngine::Create(
            spoa::AlignmentType::kNW, 3, -5, -3);
    //spoa::createAlignmentEngine(static_cast<spoa::AlignmentType>(0), 5, -10, -4, -4);

    //auto graph = spoa::createGraph();
    spoa::Graph graph{};

    for (int i = 0; i < W.size(); i++) {
        auto alignment = alignment_engine->Align(W[i], graph);
        graph.AddAlignment(alignment, W[i]);
    }

    std::vector<std::string> msa = graph.GenerateMultipleSequenceAlignment();

    return msa;
}



vector<string> easy_consensus(vector<string> V, unsigned maxMSA, string path){
    uint32_t non_empty(0);
    //~ absoluteMAJ_consensus(V);
    if(V.size()==1){
        return V;
    }
    //~ uint32_t maximum(0);
    std::set<std::string> mySet;
    for(uint32_t iV(0);iV<V.size();++iV){
        //~ cerr<<V[iV]<<endl;
        //~ maximum=max(maximum,(uint32_t)V[iV].size());
        //~ if(V[iV].size()!=0){
        //~ non_empty++;
        //~ continue;
        //~ }
        mySet.insert(V[iV]);
    }
    // if(V[iV].size()!=V[0].size()){
    if(mySet.size() > 1) {
        // std::cerr << "go consensus_POA" << std::endl;
        V =consensus_SPOA(V, maxMSA, path);

        // for (auto s : Vv) {
        //      std::cerr << s << std::endl;
        // }
        // std::cerr << std::endl;
        // for (auto s : Vvv) {
        //      std::cerr << s << std::endl;
        // }
        // std::cerr << std::endl;
        // std::cerr << std::endl;

        // return {V[0]};

        // std::cerr << "ok" << std::endl;
        // break;
    } else {
        return {V[0]};
    }
    //~ cerr<<non_empty<<"ne";
    //~ cin.get();
    string result;
    //~ for(uint i(0);i<V.size();++i){
    //~ cerr<<V[i]<<endl;
    //~ }
    //~ cerr<<"END PP"<<endl;
    for(uint32_t iS(0);iS<V[0].size();++iS){
        uint32_t cA,cC,cG,cT,cM;
        cM=cA=cC=cG=cT=0;
        for(uint32_t iV(0);iV<V.size();++iV){
            if(V[iV].size()==0){
                continue;
            }
            switch(V[iV][iS]){
                case 'A': ++cA;break;
                case 'C': ++cC;break;
                case 'G': ++cG;break;
                case 'T': ++cT;break;
                default:
                    cM++;
                    //~ cerr<<"NOPE"<<V[iV][iS]<<"?"<<endl;
                    //~ cerr<<iS<<" "<<V[iV].size()<<" "<<iV<<" "<<V.size()<<endl;
            }
        }
        if(cM>cA and cM>cC and cM>cT and cM>cG){
            // result+=('-');
            continue;
        }
        if(cA>cC and cA>cG and cA>cT){
            result+=('A');
            continue;
        }
        if(cC>cA and cC>cG and cC>cT){
            result+=('C');
            continue;
        }
        if(cG>cA and cG>cC and cG>cT){
            result+=('G');
            continue;
        }
        if(cT>cA and cT>cG and cT>cC){
            result+=('T');
            continue;
        }
        if (V[0][iS] != '-') {
            result+=(V[0][iS]);
        }
        // result+='N';
        continue;
        //~ cerr<<"TIE"<<endl;
        return V;
    }
    //~ cerr<<"EASYCONSENSU end"<<endl;

    return {result};
}

std::string concatNucR(std::string f, int i) {
    switch (i) {
        case 0:
            return f + "A";
        case 1:
            return f + "C";
        case 2:
            return f + "G";
        default:
            return f + "T";
    }
}

std::vector<std::string> getNeighbours(std::string kMer, unsigned merSize, int left, robin_hood::unordered_map<kmer, unsigned> merCounts, unsigned solidThresh) {
    std::vector<std::string> neighbours;
    std::string f, n, t = "";
    kmer k;
    std::transform(kMer.begin(), kMer.end(), kMer.begin(), ::toupper);

    if (left == 1) {
        kMer = rev_comp::run(kMer);
    }
    f = kMer.substr(1);
    int i = 0;
    for (i = 0; i < 4; i++) {
        k = str2num(f);
        k <<= 2;
        k += i;
        if (left == 1) {
            t = kmer2str(k, merSize);
            t = rev_comp::run(t);
            k = str2num(t);
        }
        if (merCounts[k] >= solidThresh) {
            if (left == 1) {
                neighbours.push_back(t);
            } else {
                neighbours.push_back(concatNucR(f, i));
            }
        }
    }

    // Sort in ascending order of number of occurrences
    std::sort(neighbours.begin(), neighbours.end(),
              [&merCounts](std::string& n1, std::string& n2) {
                  return  merCounts[str2num(n1)] > merCounts[str2num(n2)];
              }
    );
    return neighbours;
}

unsigned extendLeft(robin_hood::unordered_map<kmer, unsigned> merCounts, unsigned curK, unsigned extLen, string &LR, unsigned solidThresh) {
    vector<string> neighbours;
    vector<string>::iterator it;
    unsigned dist = 0;

    // Get the leftmost k-mer and search for a path in the graph
    neighbours = getNeighbours(LR.substr(0, curK), curK, 1, merCounts, solidThresh);
    it = neighbours.begin();

    // Keep traversing the graph while the long reads's border or a branching path aren't reached
    while (neighbours.size() == 1 && it != neighbours.end() && dist < extLen) {
        LR = (*it).substr(0, it->length() - (curK - 1)) + LR;
        dist = dist + it->length() - (curK - 1);
        // Get the leftmost k-mer and search for a path in the graph
        neighbours = getNeighbours(LR.substr(0, curK), curK, 1, merCounts, solidThresh);
        it = neighbours.begin();
    }

    return dist;
}

unsigned extendRight(robin_hood::unordered_map<kmer, unsigned> merCounts, unsigned curK, unsigned extLen, string &LR, unsigned solidThresh) {
    vector<string> neighbours;
    vector<string>::iterator it;
    unsigned dist = 0;

    // Get the leftmost k-mer and search for a path in the graph
    neighbours = getNeighbours(LR.substr(LR.length() - curK), curK, 0, merCounts, solidThresh);
    it = neighbours.begin();

    // Keep traversing the graph while the long reads's border or a branching path aren't reached
    while (it != neighbours.end() && dist < extLen) {
        LR = LR + (*it).substr(curK - 1);
        dist = dist + it->length() - (curK - 1);
        // Get the leftmost k-mer and search for a path in the graph
        neighbours = getNeighbours(LR.substr(LR.length() - curK), curK, 0, merCounts, solidThresh);
        it = neighbours.begin();
    }

    return dist;
}


int link(robin_hood::unordered_map<kmer, unsigned> merCounts, std::string srcSeed, std::string tgtSeed, unsigned curK, std::set<std::string> &visited, unsigned* curBranches, unsigned dist, std::string curExt, std::string &missingPart, unsigned merSize, unsigned LRLen, unsigned maxBranches, unsigned solidThresh, unsigned minOrder) {
    if (curK < minOrder || *curBranches > maxBranches || dist > LRLen) {
        missingPart = std::string();
        return 0;
    }

    std::string srcAnchor = curExt.substr(curExt.length() - curK);
    std::string tgtAnchor = tgtSeed.substr(0, curK);
    std::vector<std::string> neighbours;
    std::vector<std::string>::iterator it;
    bool found = srcAnchor == tgtAnchor;
    std::string curRead;
    std::string resPart1 = std::string(curExt);
    std::set<std::string>::iterator itf;

    // Search for a path in the graph starting from the source's anchor
    neighbours = getNeighbours(srcAnchor.substr(srcAnchor.length() - curK), curK, 0, merCounts, solidThresh);
    it = neighbours.begin();

    // While the destination or a braching path aren't reached, keep on traversing the graph
    while (!found && neighbours.size() == 1 && it != neighbours.end() && dist <= LRLen) {
        curRead = *it;
        itf = visited.find(curRead);
        tgtAnchor = tgtSeed.substr(0, curRead.length());
        found = curRead == tgtAnchor;
        if (!found && (itf == visited.end())) {
            visited.insert(curRead);
            resPart1 = resPart1 + curRead[curK - 1];
            dist = dist + curRead.length() - (curK - 1);

            // Update the current k-mer, and search for a path in the graph
            srcAnchor = resPart1.substr(resPart1.length() - curK);
            neighbours = getNeighbours(srcAnchor.substr(srcAnchor.length() - curK), curK, 0, merCounts, solidThresh);
            it = neighbours.begin();
        } else if (found) {
            resPart1 = resPart1 + curRead[curK - 1];
        } else {
            it++;
        }
    }

    // If a branching path is reached, explore the different possible paths with backtracking
    while (!found && neighbours.size() > 1 && it != neighbours.end() && dist <= LRLen) {
        curRead = *it;
        itf = visited.find(curRead);
        tgtAnchor = tgtSeed.substr(0, curRead.length());
        found = curRead == tgtAnchor;
        if (!found && (itf == visited.end())) {
            visited.insert(curRead);
            (*curBranches)++;
            found = link(merCounts, srcSeed, tgtSeed, merSize, visited, curBranches, dist + curRead.length() - (curK - 1), resPart1 + curRead[curK - 1], missingPart, merSize, LRLen, maxBranches, solidThresh, minOrder);
            if (!found) {
                ++it;
            } else {
                return 1;
            }
        } else if (found) {
            resPart1 = resPart1 + curRead[curK - 1];
        } else {
            ++it;
        }
    }

    // If the source couldn't be linked to the destination, try again with a graph of smaller order, otherwhise update the missing part and return
    if (!found) {
        return 0;
    } else {
        missingPart = resPart1 + tgtSeed.substr(curK);
        return 1;
    }
}


vector<vector<string>> global_consensus(const  vector<vector<string>>& V, uint32_t n, unsigned maxMSA, string path){
    vector<vector<string>> result;
    string stacked_consensus;
    for(uint32_t iV(0);iV<V.size();++iV){
        if(V[iV].size()==0){
            // cerr<<"MISSING WINDOWS"<<endl;
            continue;
        }
        // std::cerr << "go easy_consensus" << std::endl;
        vector<string> consensus(easy_consensus(V[iV], maxMSA, path));
        // std::cerr << "ok" << std::endl;
        //~ cerr<<"EASYCONSENSUS"<<endl;
        //~ cerr<<consensus[0]<<endl;
        //~ if(consensus.size()==1){
        stacked_consensus+=consensus[0];
        //~ }else{
        //~ if(stacked_consensus.size()!=0){
        //~ vector<string> vect(n, stacked_consensus);
        //~ result.push_back(vect);
        //~ stacked_consensus="";
        //~ }
        //~ result.push_back(consensus);
        //~ }
    }
    if(stacked_consensus.size()!=0){
        // stacked_consensus.erase (remove(stacked_consensus.begin(), stacked_consensus.end(), '-'), stacked_consensus.end());
        vector<string> vect(1, stacked_consensus);
        result.push_back(vect);
        stacked_consensus="";
    }
    return result;
}




std::pair<std::vector<std::vector<std::string>>, robin_hood::unordered_map<kmer, unsigned>> MSABMAAC(const vector<string>& Reads,uint32_t k, double edge_solidity, unsigned solidThresh, unsigned minAnchors, unsigned maxMSA, string path){
    int kmer_size(k);
    //~ vector<string> VTest;;
    //~ VTest.push_back("CTGACTGACCCCGTACGTCA");
    //~ VTest.push_back("CTGACTGATTTCGTACGTCA");
    //~ VTest.push_back("CTGACTGAAAACGTACGTCA");
    //~ VTest.push_back("CTGACTGAAAACGTACGTCA");
    //~ VTest.push_back("CTGACTGAAAACGTACGTCAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");
    //~ VTest.push_back("CTGACTGATTTCGTACGTCA");
    //~ VTest.push_back("CTGACTGATTTCGTACGTCA");
    //~ VTest.push_back("CTGACTGATTTCGTACGTCA");
    //~ VTest.push_back("CTGACTGATTTCGTACGTCA");
    //~ VTest.push_back("CTGACTGATTTCGTACGTCA");
    //~ VTest.push_back("CTGACTGATTTCGTACGTCA");
    //~ VTest.push_back("CTGACTGCCCCCGTACGTCA");
    //~ VTest.push_back("CTGACTGATTTCGTACGTCA");
    //~ VTest.push_back("CTGACTGATTTCGTACGTCA");
    //~ VTest.push_back("CTGACTGACCCCGTACGTCA");
    //~ auto nadine({VTest});
    //~ auto GC(global_consensus(nadine,1));
    //~ cerr<<GC[0][0]<<endl;
    //~ exit(0);



    kmer2localisation kmer_index;
    robin_hood::unordered_map<kmer, unsigned> merCounts;
    // std::cerr << "1" << std::endl;
    fill_index_kmers(Reads,kmer_index,kmer_size,merCounts, solidThresh);
    // std::cerr << "ok" << std::endl;
    // cerr<<"PHASE 1 done"<<endl;
    //~ return {};

    // std::cerr << "2" << std::endl;
    auto kmer_count(filter_index_kmers(kmer_index,edge_solidity));
    // std::cerr << "ok" << std::endl;

    // clean_suspcious_reads(kmer_index,Reads.size(),50);
    //~ auto kmer_count(filter_index_kmers(kmer_index,percent_shared));
    //~ cerr<<"PHASE 2.1 done"<<endl;
    // std::cerr << "3" << std::endl;
    auto template_read(get_template(kmer_index,Reads[0],kmer_size));
    // std::cerr << "ok" << std::endl;
    //~ cerr<<"PHASE 2 done"<<endl;

    // std::cerr << "4" << std::endl;
    vector<kmer> anchors(longest_ordered_chain(kmer_index, template_read,edge_solidity));
    // std::cerr << "ok" << std::endl;
    //~ cerr<<"PHASE 3 done"<<endl;

    // std::cerr << "5" << std::endl;
    vector<double> relative_positions=(average_distance_next_anchor(kmer_index,anchors,kmer_count,false));
    // std::cerr << "ok" << std::endl;
    //~ cerr<<"PHASE 4 done"<<endl;

    // std::cerr << "6" << std::endl;
    vector<vector<string>> result(split_reads(anchors,relative_positions,Reads,kmer_index,kmer_size));
    // std::cerr << "splits : " << result.size() << std::endl;
    if (result.size() < minAnchors) {
        // std::cerr << "to few anchors" << std::endl;
        // std::cerr << "anchors nb : " << result.size() << std::endl;
        // std::cerr << "support : " << Reads.size() << std::endl;
        std::vector<std::string> res;
        res.push_back("");
        std::vector<std::vector<std::string>> fRes;
        // fRes.push_back(res);
        return std::make_pair(fRes, merCounts);
    }
    // std::cerr << "ok" << std::endl;
    //~ cerr<<"PHASE 5 done"<<endl;


    // cerr<<""<<result.size()<<"       ";
    //~ for(uint i(0);i<result.size();++i){
    //~ for(uint j(0);j<result[i].size();++j){
    //~ cerr<<result[i][j]<<" ";
    //~ }
    //~ cerr<<endl;
    //~ }
    //~ cin.get();
    // vector<vector<string>> result;
    // result.push_back(Reads);
    // std::cerr << "7" << std::endl;
    result=global_consensus(result,Reads.size(), maxMSA, path);
    // std::cerr << "ok" << std::endl;
    //~ cerr<<"PHASE 6 done"<<endl;

    return std::make_pair(result, merCounts);
}


//
// Created by loren on 11/03/2025.
//
std::string toUpperCase(std::string& s, int beg, int end) {
    std::string res = s;
    std::locale loc;
    for (int i = beg; i <= end; i++) {
        res[i] = std::toupper(res[i], loc);
    }

    return res;
}

bool isUpperCase(char c) {
    return 'A' <= c and c <= 'Z';
}

std::string toLowerCase(std::string& s, int beg, int end) {
    std::string res = s;
    std::locale loc;
    for (int i = beg; i <= end; i++) {
        res[i] = std::tolower(res[i], loc);
    }

    return res;
}

std::string trimRead(std::string correctedRead, unsigned merSize) {
    unsigned beg, end, n;
    unsigned i;
    i = 0;
    n = 0;
    while (i < correctedRead.length() and n < merSize) {
        if (isUpperCase(correctedRead[i])) {
            n++;
        } else {
            n = 0;
        }
        i++;
    }
    beg = i - merSize;

    i = correctedRead.length() - 1;
    n = 0;
    while (i >= 0 and n < merSize) {
        if (isUpperCase(correctedRead[i])) {
            n++;
        } else {
            n = 0;
        }
        i--;
    }
    end = i + merSize;

    if (end > beg) {
        return correctedRead.substr(beg, end - beg + 1);
    } else {
        return "";
    }
}

int nbCorBases(std::string correctedRead) {
    int n = 0;
    for (unsigned i = 0; i < correctedRead.length(); i++) {
        if ('A' <= correctedRead[i] && correctedRead[i] <= 'Z') {
            n++;
        }
    }

    return n;
}

bool dropRead(std::string correctedRead) {
    return (float) nbCorBases(correctedRead) / correctedRead.length() < 0.1;
}

unsigned* getCoverages(std::string template_read, std::vector<Overlap> alignments) {
    unsigned tplLen = template_read.length();
    unsigned* coverages = (unsigned*) calloc(tplLen, sizeof(int));
    unsigned beg, end;
    unsigned i;


    return coverages;
}


std::vector<std::pair<unsigned,unsigned >> getAlignmentWindowsPositions(std::string template_read, std::vector<Overlap> alignments, unsigned minSupport,  unsigned windowSize, int overlappingWindows) {
    unsigned* coverages = getCoverages(template_read,alignments);
    unsigned i;
    unsigned beg = 0;
    unsigned tplen=template_read.length();
    unsigned end = tplen - 1;

    std::vector<std::pair<unsigned, unsigned>> pilesPos;
    unsigned curLen = 0;
    beg = 0;
    end = 0;
    i = 0;

    while (i < tplen) {
        if (curLen >= windowSize) {
            pilesPos.push_back(std::make_pair(beg, beg + curLen - 1));
            if (overlappingWindows) {
                i = i - overlappingWindows;
            }
            beg = i;
            curLen = 0;
        }
        if (coverages[i] < minSupport) {
            curLen = 0;
            i++;
            beg = i;
        } else {
            curLen++;
            i++;
        }
    }

    // Special case for the last window
    int pushed = 0;
    beg = 0;
    end = tplen - 1;
    curLen = 0;
    i = tplen - 1;
    while (i > 0 and !pushed) {
        if (curLen >= windowSize) {
            pilesPos.push_back(std::make_pair(end - curLen + 1, end));
            pushed = 1;
            end = i;
            curLen = 0;
        }
        if (coverages[i] < minSupport) {
            curLen = 0;
            i--;
            end = i;
        } else {
            curLen++;
            i--;
        }
    }

    // delete [] coverages;
    free(coverages);
    coverages = nullptr;

    return pilesPos;
}

std::vector<std::string> getAlignmentWindowsSequences(std::string template_read,std::vector<Overlap> alignments,  unsigned qBeg, unsigned end, unsigned merSize) {
    std::vector<std::string> curPile;
    std::vector<unsigned> curScore;
    unsigned length, shift;
    length = end - qBeg + 1;
    unsigned curPos = 0;
    unsigned tBeg, tEnd;

    if (qBeg + length - 1 >= template_read.length()) {
        return curPile;
    }


    curPile.push_back(template_read.substr(qBeg, length));



    std::string tmpSeq;


    for(const auto overlap:alignments){
        tBeg = overlap.target_start;
        tEnd = overlap.target_end;
        length = end - qBeg + 1;
        if (qBeg > overlap.query_start) {
            shift = qBeg - overlap.query_start;
        } else {
            shift = 0;
        }

        // For all alignments than span, or begin/end in the query window
        if ( ((overlap.query_start <= qBeg and overlap.query_end > qBeg) or (end <= overlap.query_end and overlap.query_start < end)) and overlap.target_start + shift <= overlap.target_end) {

            if (qBeg < overlap.query_start and overlap.query_end < end) {
                shift = 0;
                tBeg = std::max(0, (int) overlap.target_start - ((int) overlap.query_start - (int) qBeg));
                tEnd = std::min((int) overlap.target.length() - 1, (int) overlap.target_end + ((int) end - (int) overlap.query_end));
                length = tEnd - tBeg + 1;
            } else if (qBeg < overlap.query_start) {
                shift = 0;
                tBeg = std::max(0, (int) overlap.target_start - ((int) overlap.query_start - (int) qBeg));
                length = std::min((int) length, std::min((int) overlap.target.length() - 1, (int) tBeg + (int) length - 1) - (int) tBeg + 1);
            } else if (overlap.query_end < end) {
                tEnd = std::min((int) overlap.target.length() - 1, (int) overlap.target_end + ((int) end - (int) overlap.query_end));
                length = std::min((int) length, (int) tEnd - std::max(0, (int) tEnd - (int) length + 1) + 1);
            }


            tmpSeq = overlap.target.substr(tBeg, tEnd - tBeg + 1);

            if (overlap.strand=="-") {
                tmpSeq = rev_comp::run(tmpSeq);
            }

            tmpSeq = tmpSeq.substr(shift, length);


            // Default
            if (tmpSeq.length() >= merSize) {
                curPile.push_back(tmpSeq);
            }
        }
        curPos++;

    }

    return curPile;
}



char rev_comp::complement[int('t') + 1];
rev_comp* rev_comp::_instance = nullptr;

std::string rev_comp::run(std::string seq) {
    rev_comp::build_instance();

    auto first = seq.begin(), last = seq.end();

    while(true) {
        if(first == last || first == --last) {
            if(seq.length() % 2) {
                *first = rev_comp::complement[(unsigned char) *first];
            }
            return seq;
        } else {
            *first = rev_comp::complement[(unsigned char) *first];
            *last = rev_comp::complement[(unsigned char) *last];
            std::iter_swap(first, last);
            ++first;
        }
    }
}

void rev_comp::build_instance() {
    if(_instance == nullptr) {
        _instance = new rev_comp();
    }
}

rev_comp::rev_comp() {
    this->complement['A'] = 'T';
    this->complement['T'] = 'A';
    this->complement['C'] = 'G';
    this->complement['G'] = 'C';
    this->complement['a'] = 't';
    this->complement['t'] = 'a';
    this->complement['c'] = 'g';
    this->complement['g'] = 'c';
}

rev_comp::~rev_comp(){

}

robin_hood::unordered_map<std::string, std::vector<unsigned>> getKMersPos(std::string sequence, unsigned merSize) {
    robin_hood::unordered_map<std::string, std::vector<unsigned>> mers;

    for (unsigned i = 0; i < sequence.length() - merSize + 1; i++) {
        mers[sequence.substr(i, merSize)].push_back(i);
    }

    return mers;
}

int getNextSrc(std::string correctedRead, unsigned beg, unsigned merSize) {
    unsigned nb = 0;
    unsigned i = beg;

    while (i < correctedRead.length() and (isUpperCase(correctedRead[i]) or nb < merSize)) {
        if (isUpperCase(correctedRead[i])) {
            nb++;
        } else {
            nb = 0;
        }
        i++;
    }

    return nb >= merSize ? i - 1 : -1;
}

int getNextDst(std::string correctedRead, unsigned beg, unsigned merSize) {
    unsigned nb = 0;
    unsigned i = beg;

    while (i < correctedRead.length() and nb < merSize) {
        if (isUpperCase(correctedRead[i])) {
            nb++;
        } else {
            nb = 0;
        }
        i++;
    }

    return nb >= merSize ? i - 1 : -1;
}


// Anchors without repeated k-mers
std::vector<std::pair<std::string, std::string>> getAnchors(robin_hood::unordered_map<kmer, unsigned>& merCounts, std::string srcZone, std::string dstZone, unsigned merSize, unsigned nb) {
    std::vector<std::pair<std::string, std::string>> res;
    unsigned i;

    robin_hood::unordered_map<std::string, std::vector<unsigned>> mersPosSrc = getKMersPos(srcZone, merSize);
    robin_hood::unordered_map<std::string, std::vector<unsigned>> mersPosDst = getKMersPos(dstZone, merSize);

    // Consider all k-mers of the src zone as potential anchors
    std::vector<std::string> candidatesSrc(srcZone.size() - merSize + 1);
    for (i = 0; i < srcZone.size() - merSize + 1; i++) {
        candidatesSrc[i] = srcZone.substr(i, merSize);
    }
    // Same with the dst zone
    std::vector<std::string> candidatesDst(dstZone.size() - merSize + 1);
    for (i = 0; i < dstZone.size() - merSize + 1; i++) {
        candidatesDst[i] = dstZone.substr(i, merSize);
    }

    // Add the anchors pairs to the result vector, without allowing repeated k-mers
    for (std::string csrc : candidatesSrc) {
        if (mersPosSrc[csrc].size() == 1) {
            for (std::string cdst : candidatesDst) {
                if (mersPosDst[cdst].size() == 1) {
                    res.push_back(std::make_pair(csrc, cdst));
                }
            }
        }
    }

    // Sort the anchors vector in ascending order of the number of occurrences of (src + dst)
    std::sort(res.begin(), res.end(),
              [&merCounts](std::pair<std::string, std::string>& r1, std::pair<std::string, std::string>& r2) {
                  int occ1 = merCounts[str2num(r1.first)] + merCounts[str2num(r1.second)];
                  int occ2 = merCounts[str2num(r2.first)] + merCounts[str2num(r2.second)];
                  return occ1 > occ2;
              }
    );

    std::vector<std::pair<std::string, std::string>> finalRes;
    for (i = 0; i < nb and i < res.size(); i++) {
        finalRes.push_back(res[i]);
    }

    return finalRes;
}

std::string polishCorrection(std::string correctedRead, robin_hood::unordered_map<kmer, unsigned>& merCounts, unsigned merSize, int solidThresh) {
    std::set<std::string> visited;
    unsigned curBranches;
    unsigned dist;
    std::string curExt;
    std::string correctedRegion;
    unsigned maxSize;
    unsigned maxBranches = 50;
    std::vector<std::pair<std::string, std::string>> corList;
    int zone = 3;
    int srcBeg, srcEnd, dstBeg, dstEnd;
    unsigned tmpSrcBeg = 0, tmpSrcEnd = 0, tmpDstBeg = 0, tmpDstEnd = 0;
    std::string src, dst;
    std::pair<int, int> pos;
    std::vector<std::pair<std::string, std::string>> anchors;
    unsigned anchorNb;
    std::string srcZone, dstZone;
    robin_hood::unordered_map<std::string, std::vector<unsigned>> srcPos, dstPos;
    std::string oldCorrectedRead;
    int b, l;
    std::string r, c;

    // Skip uncorrected head of the read
    unsigned i = 0;
    while (i < correctedRead.length() and !isUpperCase(correctedRead[i])) {
        i++;
    }

    if (i > 0 and i < correctedRead.length() and correctedRead.length() - i >= merSize) {
        int extLen = i;
        oldCorrectedRead = correctedRead;
        correctedRead = correctedRead.substr(i);
        int extSize = extendLeft(merCounts, merSize, extLen, correctedRead, solidThresh);
        if (extSize < extLen) {
            correctedRead = oldCorrectedRead.substr(0, extLen - extSize) + correctedRead;
            i = i - (extLen - extSize);
        }
    }

    // Search for poorly supported regions bordered by solid corrected regions
    while (i < correctedRead.length()) {
        srcEnd = getNextSrc(correctedRead, i, merSize + zone);
        dstEnd = getNextDst(correctedRead, srcEnd + 1, merSize + zone);
        srcBeg = srcEnd - merSize - zone + 1;
        dstBeg = dstEnd - merSize - zone + 1;

        // Polish the poorly supported region region if 2 anchors were found
        if (srcEnd != -1 and dstEnd != -1) {
            correctedRegion = "";
            srcZone = correctedRead.substr(srcBeg, merSize + zone);
            dstZone = correctedRead.substr(dstBeg, merSize + zone);
            anchors = getAnchors(merCounts, srcZone, dstZone, merSize, 5);
            srcPos = getKMersPos(srcZone, merSize);
            dstPos = getKMersPos(dstZone, merSize);

            // Attempt to link frequent anchors
            anchorNb = 0;
            while (anchorNb < anchors.size() and correctedRegion.empty()) {
                src = anchors[anchorNb].first;
                dst = anchors[anchorNb].second;
                tmpSrcBeg = srcBeg + srcPos[src][0];
                tmpSrcEnd = tmpSrcBeg + merSize - 1;
                tmpDstBeg = dstBeg + dstPos[dst][0];
                tmpDstEnd = tmpDstBeg + merSize - 1;

                if (src != dst) {
                    curBranches = 0;
                    dist = 0;
                    curExt = src;
                    correctedRegion = "";
                    maxSize = 15.0 / 100.0 * 2.0 * (tmpDstBeg - tmpSrcEnd - 1) + (tmpDstBeg - tmpSrcEnd - 1) + merSize;
                    link(merCounts, src, dst, merSize, visited, &curBranches, dist, curExt, correctedRegion, merSize, maxSize, maxBranches, solidThresh, merSize);
                }
                anchorNb++;
            }

            if (!correctedRegion.empty()) {
                // Anchor the correction to the read
                r = correctedRead.substr(tmpSrcBeg, tmpDstEnd - tmpSrcBeg + 1);
                c = correctedRegion;
                b = correctedRead.find(r);
                l = r.length();
                if ((int) b != -1) {
                    correctedRead.replace(b, l, c);
                    i = b;
                } else {
                    i = tmpDstBeg > i ? tmpDstBeg : dstBeg;
                }
            } else {
                i = tmpDstBeg > i ? tmpDstBeg : dstBeg;
            }
        } else {
            i = correctedRead.length();
        }
    }

    i = correctedRead.length() - 1;
    while (i > 0 and !isUpperCase(correctedRead[i])) {
        i--;
    }

    if (i > 0 and i < correctedRead.length() - 1 and i + 1 >= merSize) {
        int extLen = correctedRead.length() - 1 - i;
        oldCorrectedRead = correctedRead;
        correctedRead = correctedRead.substr(0, i + 1);
        int extSize = extendRight(merCounts, merSize, extLen, correctedRead, solidThresh);
        if (extSize < extLen) {
            correctedRead = correctedRead + oldCorrectedRead.substr(oldCorrectedRead.length() - (extLen - extSize), extLen - extSize);
        }
    }

    return correctedRead;
}

int nbSolidMers(std::string seq, robin_hood::unordered_map<kmer, unsigned> merCounts, unsigned merSize, unsigned solidThresh) {
    int nb = 0;
    for (unsigned i = 0; i < seq.length() - merSize + 1; i++) {
        if (merCounts[str2num(seq.substr(i, merSize))] >= solidThresh) {
            nb++;
        }
    }

    return nb;
}

int nbUpperCase(std::string s) {
    int nb = 0;
    for (unsigned i = 0; i < s.length(); i++) {
        if (isUpperCase(s[i])) {
            nb++;
        }
    }

    return nb;
}

std::pair<int, int> getIndels(std::string cigar){
    int ins = 0;
    int del = 0;
    int current = 0;
    for(unsigned i = 0; i < cigar.length(); i++){
        if('0' <= cigar[i] && cigar[i] <= '9'){
            current = (current * 10) + (cigar[i] - '0');
        } else {
            if (cigar[i] == 'I') {
                ins += current;
            } else if (cigar[i] == 'D') {
                del += current;
            }
            current = 0;
        }
    }
    return std::make_pair(ins, del);

}

std::string alignConsensus(std::string sequence, std::vector<std::string>& consensuses, std::vector<robin_hood::unordered_map<kmer, unsigned>>& merCounts, std::vector<std::pair<unsigned, unsigned>>& pilesPos, std::vector<std::string>& templates, int startPos, unsigned windowSize, unsigned windowOverlap, unsigned solidThresh, unsigned merSize) {
    StripedSmithWaterman::Aligner aligner;
    StripedSmithWaterman::Filter filter;
    StripedSmithWaterman::Alignment alignment;
    StripedSmithWaterman::Alignment subAlignment;
    int32_t maskLen = 15;
    unsigned beg, end, oldEnd;
    oldEnd = 0;
    std::string outSequence;
    outSequence = sequence;
    std::transform(outSequence.begin(), outSequence.end(), outSequence.begin(), ::tolower);

    std::string corWindow;
    unsigned i = 0;
    std::string tmpSequence, consUp;
    int curPos = startPos;
    int alPos;
    int sizeAl;
    std::string curCons, oldCons;
    robin_hood::unordered_map<kmer, unsigned> oldMers;
    robin_hood::unordered_map<kmer, unsigned> curMers;
    unsigned overlap;
    std::string seq1, seq2;
    int solidMersSeq1, solidMersSeq2;
    std::pair<int, int> indels;
    unsigned ins, del;

    for (i = 0; i < consensuses.size(); i++) {
        if (consensuses[i].length() < merSize) {
            curCons = templates[i];
            std::transform(consUp.begin(), consUp.end(), consUp.begin(), ::tolower);
        } else {
            curCons = consensuses[i];
        }
        curMers = merCounts[i];

        alPos = std::max(0, (int) curPos - (int) windowOverlap);
        if (alPos + windowSize + 2 * windowOverlap >= outSequence.length()) {
            sizeAl = outSequence.length() - alPos;
        } else {
            sizeAl = windowSize + 2 * windowOverlap;
        }

        aligner.Align(curCons.c_str(), outSequence.c_str() + alPos, sizeAl, filter, &alignment, maskLen);
        beg = alignment.ref_begin + alPos;
        end = alignment.ref_end + alPos;
        curCons = curCons.substr(alignment.query_begin, alignment.query_end - alignment.query_begin + 1);

        // Check if alignment positions overlap the previous window. If they do, chose the best subsequence
        if (i != 0 and oldEnd >= beg) {
            overlap = oldEnd - beg + 1;
            if (consensuses[i].length() >= merSize and oldCons.length() >= overlap and curCons.length() >= overlap) {
                seq1 = oldCons.substr(oldCons.length() - 1 - overlap + 1, overlap);
                seq2 = curCons.substr(0, overlap);
                if (toUpperCase(seq1, 0, seq1.length() - 1) != toUpperCase(seq2, 0, seq2.length() - 1)) {
                    if (overlap >= merSize) {
                        solidMersSeq1 = nbSolidMers(seq1, oldMers, merSize, solidThresh);
                        solidMersSeq2 = nbSolidMers(seq2, curMers, merSize, solidThresh);
                    } else {
                        solidMersSeq1 = nbUpperCase(seq1);
                        solidMersSeq2 = nbUpperCase(seq2);
                    }
                    if (solidMersSeq1 > solidMersSeq2) {
                        aligner.Align(seq1.c_str(), seq2.c_str(), std::min(seq1.length(), seq2.length()), filter, &subAlignment, maskLen);
                        indels = getIndels(subAlignment.cigar_string);
                        ins = indels.first;
                        del = indels.second;
                        if (overlap - ins + del < curCons.length()) {
                            curCons = seq1 + curCons.substr(overlap - ins + del);
                        } else {
                            curCons = "";
                        }
                    }
                }
            }
        }

        if (curCons != "") {
            if (consensuses[i].length() >= merSize){
                consUp = curCons;
                std::transform(consUp.begin(), consUp.end(), consUp.begin(), ::toupper);
                outSequence.replace(beg, end - beg + 1, consUp);
            }
            if (i < consensuses.size() - 1) {
                curPos = curPos + pilesPos[i+1].first - pilesPos[i].first - (end - beg + 1) + curCons.length() ;
                oldCons = curCons;
                oldMers = merCounts[i];
                oldEnd = beg + curCons.length() - 1;
            }
        }
    }

    return outSequence;
}


std::string weightConsensus(std::string& consensus, std::vector<std::string>& pile, robin_hood::unordered_map<kmer, unsigned>& merCounts, unsigned merSize, unsigned windowSize, unsigned solidThresh) {
    std::vector<std::string> splits;
    std::string curSplit;

    std::string header = "";
    std::string sequence = "";
    std::string curFct;

    unsigned i = 0;
    while (i < consensus.length() - merSize + 1) {
        curFct = consensus.substr(i, merSize);
        curFct = toUpperCase(curFct, 0, merSize);
        if (merCounts[str2num(curFct)] >= solidThresh) {
            consensus = toUpperCase(consensus, i, i + merSize - 1);
        } else {
            consensus = toLowerCase(consensus, i, i + merSize - 1);
        }
        i++;
    }

    return consensus;
}

std::pair<std::string, robin_hood::unordered_map<kmer, unsigned>> computeConsensusReadCorrection(std::vector<std::string> & piles, std::pair<unsigned, unsigned>& pilesPos, unsigned& minSupport, unsigned& merSize, unsigned& commonKMers, unsigned& minAnchors, unsigned& solidThresh, unsigned& windowSize, unsigned maxMSA, std::string path) {
    int bmeanSup;
    bmeanSup = std::min((int) commonKMers, (int) piles.size() / 2);
    std::pair<std::vector<std::vector<std::string>>, robin_hood::unordered_map<kmer, unsigned>> rOut = MSABMAAC(piles, merSize, bmeanSup, solidThresh, minAnchors, maxMSA, path);

    if (rOut.first.size() == 0) {
        return std::make_pair(piles[0], rOut.second);
    }
    auto result = rOut.first;
    auto merCounts = rOut.second;
    std::string corTpl = result[0][0];

    // Polish the consensus
    std::vector<std::pair<std::string, std::string>> corList;
    if (corTpl.length() >= merSize) {
        corTpl = weightConsensus(corTpl, piles, merCounts, merSize, windowSize, solidThresh);
        corTpl = polishCorrection(corTpl, merCounts, merSize, solidThresh);
    }

    return std::make_pair(corTpl, merCounts);
}

std::string consent_correction(std::string template_read, std::vector<Overlap> alignments) {


    unsigned minSupport=0;
    unsigned minAnchors=2;
    unsigned solidThresh=2;
    unsigned windowSize=500;
    unsigned commonKMers=8;
    unsigned maxMSA=158;
    unsigned windowOverlap=50;
    bool doTrimRead=true;
    std::string readId="current_read";





    std::vector<std::pair<unsigned,unsigned >> pilesPos=getAlignmentWindowsPositions(template_read, alignments,minSupport, windowSize, windowOverlap);



    //std::cout << "Finestre di allineamento:" << std::endl;
    //for (const auto& window : pilesPos) {
    //    std::cout << "[" << window.first << ", " << window.second << "]" << std::endl;
    //}

    unsigned i = 0;
    unsigned merSize=9;

    //std::pair<std::string , robin_hood::unordered_map<kmer, unsigned>> resCons;
    //std::vector<std::string> consensuses(pilesPos.size());
    //std::vector<robin_hood::unordered_map<kmer, unsigned>> merCounts(pilesPos.size());
    //std::vector<std::string> curPile;
    //std::vector<std::string> templates(pilesPos.size());

    //for (i = 0; i < pilesPos.size(); i++) {
    //    curPile = getAlignmentWindowsSequences(template_read,alignments, pilesPos[i].first, pilesPos[i].second, 9);
    //    templates[i] = curPile[0];
    //    resCons = computeConsensusReadCorrection(curPile, pilesPos[i], minSupport, merSize, commonKMers, minAnchors, solidThresh, windowSize, maxMSA,"tmp");
    //    if (resCons.first.length() < merSize) {
    //        consensuses[i] = resCons.first;
    //    } else {
    //        consensuses[i] = resCons.first;
    //    }
    //    merCounts[i] = resCons.second;
    //}

    //std::string correctedRead = alignConsensus(template_read, consensuses, merCounts, pilesPos, templates, pilesPos[0].first, windowSize, windowOverlap, solidThresh, merSize);

    //if (doTrimRead) {
    //    correctedRead = trimRead(correctedRead, 1);
    //    if (!dropRead(correctedRead)) {
    //        std::cout<<correctedRead<<std::endl;
    //    } else {
    //        std::cout<<""<<std::endl;
    //    }
    //} else {
    //    std::cout<<correctedRead<<std::endl;
    //}

    return "";
}
