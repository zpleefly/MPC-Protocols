#include <iostream>
#include <vector>
#include <NTL/ZZ.h>
#include "OTExtension.h"

using namespace std;
using namespace NTL;

typedef OT_Pairs Labels;

enum GateType{AND, XOR};
enum PlayerType{Garbler, Evaluator};

union ID{
    int gate_id;
    struct UserID{
        int id, index;
    } user_id;
    int output_index;
};

// needs to be changed
ZZ H(ZZ k_a, ZZ k_b, int j){
    ZZ y, seed;
    seed = k_a + k_b + j;
    RandomStream initial_state = GetCurrentRandomStream();
    SetSeed(seed);
    y = RandomLen_ZZ(KAPPA);
    SetSeed(initial_state);
    return y;
}

class Wire{
public:
    bool circuit_input_wire;
    ID id;
    ZZ label;
    
    Wire(bool circuit_input_wire, bool circuit_output_wire, int id, int index = 0){
        this->circuit_input_wire = circuit_input_wire;
        if(circuit_output_wire)
            this->id.output_index = id;
        else if(circuit_input_wire){
            this->id.user_id.id = id;
            this->id.user_id.index = index;
        }
        else
            this->id.gate_id = id;
    }
};

/*
class WireLabels{
public:
    Labels k;
    bool *p;
    
    WireLabels(){};
    
    void SetLength(int n){
        this->k.x_0.SetLength(n);
        this->k.x_1.SetLength(n);
        this->p = new bool[n];
    }
};
*/

struct WireLabels{
    ZZ k_0, k_1;
    bool p;
};

struct GateLabels{
    WireLabels input_1, input_2, output;
};

/*
class GateLabels{
public:
    WireLabels input_1, input_2, output;
    
    GateLabels(int n){
        this->input_1.SetLength(n);
        this->input_2.SetLength(n);
        this->output.SetLength(n);
    }
};
*/

struct GarbledTable{
    ZZ e[2][2];
};

struct GarbledOutputTable{
    ZZ e[2];
};

class Gate{
public:
    GateType type;
    int id;
    bool output_gate;
    Wire* input_1;
    Wire* input_2;
    Wire* output;
    GarbledTable *table;
    GarbledOutputTable *output_table;
    
    Gate(){};
    
    Gate(GateType type, int id, bool output_gate = false, Wire* input_1 = NULL, Wire* input_2 = NULL, Wire* output = NULL){
        this->type = type;
        this->input_1 = input_1;
        this->input_2 = input_2;
        this->output = output;
        this->id = id;
        this->output_gate = output_gate;
        if(type == AND)
            this->table = new GarbledTable;
        if(output_gate)
            this->output_table = new GarbledOutputTable;
    }
    
    void initialize(GateType type, int id, bool output_gate = false, Wire* input_1 = NULL, Wire* input_2 = NULL, Wire* output = NULL){
        this->type = type;
        this->input_1 = input_1;
        this->input_2 = input_2;
        this->output = output;
        this->id = id;
        this->output_gate = output_gate;
        if(type == AND)
            this->table = new GarbledTable;
        if(output_gate)
            this->output_table = new GarbledOutputTable;
    }
    
    Gate(GateType type){
        this->type = type;
    }
    
    void set_inputs(Wire* input_1, Wire* input_2){
        this->input_1 = input_1;
        this->input_2 = input_2;
    }
    
    void set_output(Wire* output){
        this->output = output;
    }
};

class Player{
private:
    vector<bool> input;
    ZZ r;
    GateLabels *labels;
    Gate* gates;
    friend class GarbledCircuit;
public:
    PlayerType type;
    Vec<ZZ> x;
    OTSender* S;
    OTReceiver* R;
    
    Player(PlayerType type, vector<bool> input){
        this->type = type;
        this->input = input;
    }
};

class GarbledCircuit{
public:
    vector<bool> output;
    int n;
    Player *P_1, *P_2;
    OTExtension *ot;
    
    GarbledCircuit(Gate* gates, int number_of_gates, Player *P_1, Player *P_2){
        this->n = number_of_gates;
        this->P_1 = P_1;
        this->P_2 = P_2;
        if(this->P_1->type == Garbler)
            this->P_1->labels = new GateLabels[this->n];
        if(this->P_2->type == Evaluator){
            this->P_2->gates = gates;
            this->P_2->x.SetLength(this->P_1->input.size());
        }
        this->ot = new OTExtension((int) this->P_2->input.size(), KAPPA);
    }
    
    void garbling_phase(){
        P_1->r = RandomLen_ZZ(KAPPA - 1);
        Labels y_labels;
        int output_size = 0;
        y_labels.x_0.SetLength(P_2->input.size());
        y_labels.x_1.SetLength(P_2->input.size());
        Wire *input_1, *input_2, *output;
        for(int i = 0; i < this->n; i++){
            input_1 = P_2->gates[i].input_1;
            input_2 = P_2->gates[i].input_2;
            output = P_2->gates[i].output;
            
            if(input_1->circuit_input_wire){
                P_1->labels[i].input_1.k_0 = RandomLen_ZZ(KAPPA - 1);
                P_1->labels[i].input_1.p = RandomBits_long(1);
                P_1->labels[i].input_1.k_1 = P_1->labels[i].input_1.k_0 ^ P_1->r;
                if(input_1->id.user_id.id == 2){
                    y_labels.x_0[input_1->id.user_id.index] = (P_1->labels[i].input_1.k_0 << 1) + P_1->labels[i].input_1.p;
                    y_labels.x_1[input_1->id.user_id.index] = (P_1->labels[i].input_1.k_1 << 1) + (!P_1->labels[i].input_1.p);
                }
                else{
                    if(P_1->input[input_1->id.user_id.index])
                        P_2->x[input_1->id.user_id.index] = (P_1->labels[i].input_1.k_1 << 1) + (!P_1->labels[i].input_1.p);
                    else
                        P_2->x[input_1->id.user_id.index] = (P_1->labels[i].input_1.k_0 << 1) + P_1->labels[i].input_1.p;
                }
            }
            else{
                P_1->labels[i].input_1.k_0 = P_1->labels[input_1->id.gate_id].output.k_0;
                P_1->labels[i].input_1.k_1 = P_1->labels[input_1->id.gate_id].output.k_1;
                P_1->labels[i].input_1.p = P_1->labels[input_1->id.gate_id].output.p;
            }
            
            if(input_2->circuit_input_wire){
                P_1->labels[i].input_2.k_0 = RandomLen_ZZ(KAPPA - 1);
                P_1->labels[i].input_2.p = RandomBits_long(1);
                P_1->labels[i].input_2.k_1 = P_1->labels[i].input_2.k_0 ^ P_1->r;
                if(input_2->id.user_id.id == 2){
                    y_labels.x_0[input_2->id.user_id.index] = (P_1->labels[i].input_2.k_0 << 1) + P_1->labels[i].input_2.p;
                    y_labels.x_1[input_2->id.user_id.index] = (P_1->labels[i].input_2.k_1 << 1) + (!P_1->labels[i].input_2.p);
                }
                else{
                    if(P_1->input[input_2->id.user_id.index])
                        P_2->x[input_2->id.user_id.index] = (P_1->labels[i].input_2.k_1 << 1) + (!P_1->labels[i].input_2.p);
                    else
                        P_2->x[input_2->id.user_id.index] = (P_1->labels[i].input_2.k_0 << 1) + P_1->labels[i].input_2.p;
                }
            }
            else{
                P_1->labels[i].input_2.k_0 = P_1->labels[input_2->id.gate_id].output.k_0;
                P_1->labels[i].input_2.k_1 = P_1->labels[input_2->id.gate_id].output.k_1;
                P_1->labels[i].input_2.p = P_1->labels[input_2->id.gate_id].output.p;
            }
            
            if(P_2->gates[i].type == XOR){
                P_1->labels[i].output.k_0 = P_1->labels[i].input_1.k_0 ^ P_1->labels[i].input_2.k_0;
                P_1->labels[i].output.k_1 = P_1->labels[i].input_1.k_0 ^ P_1->labels[i].input_2.k_0 ^ P_1->r;
                P_1->labels[i].output.p = P_1->labels[i].input_1.p ^ P_1->labels[i].input_2.p;
            }
            else if(P_2->gates[i].type == AND){
                P_1->labels[i].output.k_0 = RandomLen_ZZ(KAPPA - 1);
                P_1->labels[i].output.p = RandomBits_long(1);
                P_1->labels[i].output.k_1 = P_1->labels[i].output.k_0 ^ P_1->r;
                
                // setup garbled table
                P_2->gates[i].table->e[P_1->labels[i].input_1.p][P_1->labels[i].input_2.p] = H(P_1->labels[i].input_1.k_0, P_1->labels[i].input_2.k_0, i) ^ (P_1->labels[i].output.k_0 << 1) + (P_1->labels[i].output.p);
                P_2->gates[i].table->e[P_1->labels[i].input_1.p][!P_1->labels[i].input_2.p] = H(P_1->labels[i].input_1.k_0, P_1->labels[i].input_2.k_1, i) ^ (P_1->labels[i].output.k_0 << 1) + (P_1->labels[i].output.p);
                P_2->gates[i].table->e[!P_1->labels[i].input_1.p][P_1->labels[i].input_2.p] = H(P_1->labels[i].input_1.k_1, P_1->labels[i].input_2.k_0, i) ^ (P_1->labels[i].output.k_0 << 1) + (P_1->labels[i].output.p);
                P_2->gates[i].table->e[!P_1->labels[i].input_1.p][!P_1->labels[i].input_2.p] = H(P_1->labels[i].input_1.k_1, P_1->labels[i].input_2.k_1, i) ^ (P_1->labels[i].output.k_1 << 1) + (!P_1->labels[i].output.p);
            }
            
            if(P_2->gates[i].output_gate){
                // setup garbled output table
                P_2->gates[i].output_table->e[P_1->labels[i].output.p] = H(P_1->labels[i].output.k_0, conv<ZZ>(0), i) ^ 0;
                P_2->gates[i].output_table->e[!P_1->labels[i].output.p] = H(P_1->labels[i].output.k_1, conv<ZZ>(0), i) ^ 1;
                output_size++;
            }
        }
        P_1->S = new OTSender(true, y_labels, (int) P_2->input.size(), KAPPA);
        bool* selection_bits = new bool[P_2->input.size()];
        copy(P_2->input.begin(), P_2->input.end(), selection_bits);
        P_2->R = new OTReceiver(true, selection_bits, (int) P_2->input.size());
        this->output.resize(output_size);
    }
    
    void OT_phase(){
        ot->set_participants(P_1->S, P_2->R);
        ot->base_OT();
        ot->extension_phase();
    }
    
    void evaluation_phase(){
        Wire *input_1, *input_2, *output;
        Vec<ZZ> y = P_2->R->get_selected_data();
        ZZ label, k_a, k_b, k;
        bool p_a, p_b, p;
        for(int i = 0; i < this->n; i++){
            input_1 = P_2->gates[i].input_1;
            input_2 = P_2->gates[i].input_2;
            output = P_2->gates[i].output;
            
            if(input_1->circuit_input_wire){
                if(input_1->id.user_id.id == 1)
                    input_1->label = P_2->x[input_1->id.user_id.index];
                else
                    input_1->label = y[input_1->id.user_id.index];
            }
            else
                input_1->label = P_2->gates[input_1->id.gate_id].output->label;
            
            if(input_2->circuit_input_wire){
                if(input_2->id.user_id.id == 1)
                    input_2->label = P_2->x[input_2->id.user_id.index];
                else
                    input_2->label = y[input_2->id.user_id.index];
            }
            else
                input_2->label = P_2->gates[input_2->id.gate_id].output->label;
            
            if(P_2->gates[i].type == XOR)
                output->label = input_1->label ^ input_2->label;
            else if(P_2->gates[i].type == AND){
                p_a = bit(input_1->label, 0);
                p_b = bit(input_2->label, 0);
                //cout << "p_a: " << p_a << " ; p_b: " << p_b << endl;
                k_a = input_1->label >> 1;
                k_b = input_2->label >> 1;
                output->label = H(k_a, k_b, i) ^ P_2->gates[i].table->e[p_a][p_b];
            }
            
            if(P_2->gates[i].output_gate){
                p = bit(output->label, 0);
                k = output->label >> 1;
                this->output[output->id.output_index] = conv<int>(H(k, conv<ZZ>(0), i) ^ P_2->gates[i].output_table->e[p]);
            }
        }
    }
    
    void show_views(){
        cout << "R: " << P_1->r << endl;
        for(int i = 0; i < this->n; i++){
            cout << "Gate " << i << ": " << endl;
            cout << "Input 1: k_0: " << (P_1->labels[i].input_1.k_0 << 1) + P_1->labels[i].input_1.p << "; k_1: " << (P_1->labels[i].input_1.k_1 << 1) + (!P_1->labels[i].input_1.p) << "; p: " << P_1->labels[i].input_1.p << endl;
            cout << "Input 1: k_v: " << P_2->gates[i].input_1->label << endl;
            cout << "Input 2: k_0: " << (P_1->labels[i].input_2.k_0 << 1) + P_1->labels[i].input_2.p << "; k_1: " << (P_1->labels[i].input_2.k_1 << 1) + (!P_1->labels[i].input_2.p) << "; p: " << P_1->labels[i].input_2.p << endl;
            cout << "Input 2: k_v: " << P_2->gates[i].input_2->label << endl;
            cout << "Output: k_0: " << (P_1->labels[i].output.k_0 << 1) + P_1->labels[i].output.p << "; k_1: " << (P_1->labels[i].output.k_1 << 1) + (!P_1->labels[i].output.p) << "; p: " << P_1->labels[i].output.p << endl;
            cout << "Output: k_v: " << P_2->gates[i].output->label << endl;
        }
    }
};

int main(){
    Gate* gates = new Gate[2];
    gates[0].initialize(AND, 0, true);
    gates[1].initialize(XOR, 1, true);
    Wire *input_01, *input_11, *input_02, *input_12, *output_0, *output_1;
    input_01 = new Wire(true, false, 1, 0);
    input_02 = new Wire(true, false, 2, 0);
    output_0 = new Wire(false, true, 0);
    input_12 = new Wire(false, false, 0);
    input_11 = new Wire(true, false, 1, 1);
    output_1 = new Wire(false, true, 1);
    gates[0].set_inputs(input_01, input_02);
    gates[0].set_output(output_0);
    gates[1].set_inputs(input_11, input_12);
    gates[1].set_output(output_1);
    vector<bool> x;
    vector<bool> y;
    x.push_back(true);
    y.push_back(true);
    x.push_back(true);
    Player* P_1 = new Player(Garbler, x);
    Player* P_2 = new Player(Evaluator, y);
    GarbledCircuit GC(gates, 2, P_1, P_2);
    GC.garbling_phase();
    GC.OT_phase();
    GC.evaluation_phase();
    vector<bool>::iterator it;
    int i = 0;
    cout << "Output: " << endl;
    for(it = GC.output.begin(); it != GC.output.end(); it++)
        cout << "Gate " << i++ << ": " << (*it) << endl;
    //GC.show_views();
    return 0;
}
