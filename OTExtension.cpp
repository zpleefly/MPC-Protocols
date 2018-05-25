#include <iostream>
#include <NTL/ZZ.h>
#include <NTL/vector.h>

#define SEC_LEVEL 1024
#define KAPPA 80
#define N 100000//10000000//160
#define L 80

double start, stop;

using namespace std;
using namespace NTL;

class Sender;
class Receiver;

struct OT_Pairs{
    Vec<ZZ> x_0, x_1;
};

void find_generator(ZZ& g, const ZZ& p, const ZZ& q){
    ZZ x, exp;
    sub(exp, p, 1);
    div(exp, exp, q);
    do{
        x = RandomBnd(p);
        PowerMod(g, x, exp, p);
    }while(g == 1);
    return;
}

Vec<ZZ> transpose_vector(const Vec<ZZ>& x, int n, int k){
    Vec<ZZ> y;
    y.SetLength(n);
    for(int i = 0; i < n; i++)
        for(int j = 0; j < k; j++)
            if(bit(x[j], i))
                SetBit(y[i], j);
    /*
     for(int i = 0; i < n; i++){
     for(int j = 0; j < k; j++)
     cout << bit(y[i], k - j - 1) << " ";
     cout << y[i] << endl;
     }
     cout << endl;
     for(int i = 0; i < k; i++){
     for(int j = 0; j < n; j++)
     cout << bit(x[i], n - j - 1) << " ";
     cout << x[i] << endl;
     }
     cout << endl;
     //*/
    return y;
}

ZZ ZZ_from_bits(bool* bit_sequence, int offset, int n){
    ZZ y = conv<ZZ>(0);
    for(int i = offset; i < offset + n; i++)
        if(bit_sequence[i])
            SetBit(y, i - offset);
    return y;
}

ZZ KDF(const ZZ& seed, int l){
    ZZ y;
    RandomStream initial_state = GetCurrentRandomStream();
    SetSeed(seed);
    y = RandomBits_ZZ(l);
    SetSeed(initial_state);
    return y;
}

ZZ PRG(const ZZ& seed, int l){
    ZZ y;
    RandomStream initial_state = GetCurrentRandomStream();
    SetSeed(seed);
    y = RandomBits_ZZ(l);
    SetSeed(initial_state);
    return y;
}

ZZ CRF(const ZZ& seed, int l){
    ZZ y;
    RandomStream initial_state = GetCurrentRandomStream();
    SetSeed(seed);
    y = RandomBits_ZZ(l);
    SetSeed(initial_state);
    return y;
}

class Sender{
private:
    OT_Pairs x;
public:
    int n, l;
    OT_Pairs h;
    Vec<ZZ> u, q;
    Receiver* R;
    friend class Receiver;
    friend class OT;
    friend class OTExtension;
    
    Sender(bool extended, OT_Pairs x, int n = N, int l = L){
        this->x = x;
        this->n = n;
        this->l = l;
        if(!extended){
            this->h.x_0.SetLength(n);
            this->h.x_1.SetLength(n);
        }
        else{
            this->u.SetLength(KAPPA);
            this->q.SetLength(KAPPA);
        }
    }
    
    OT_Pairs get_data(){
        return x;
    }
};

class Receiver{
private:
    bool* selection_bits;
    Vec<ZZ> x, a, t;
public:
    int n;
    OT_Pairs v;
    ZZ u;
    Sender* S;
    friend class Sender;
    friend class OT;
    friend class OTExtension;
    
    Receiver(bool extended, bool* selection_bits, int n = N){
        this->selection_bits = selection_bits;
        this->n = n;
        this->x.SetLength(n);
        if(!extended){
            this->a.SetLength(n);
            this->v.x_0.SetLength(n);
            this->v.x_1.SetLength(n);
        }
        else{
            this->t.SetLength(KAPPA);
            this->v.x_0.SetLength(KAPPA);
            this->v.x_1.SetLength(KAPPA);
        }
    }
    
    bool* get_selection_bits(){
        return selection_bits;
    }
    
    Vec<ZZ> get_selected_data(){
        return x;
    }
};

class OT{
public:
    int n, l;
    ZZ p, q, g;
    Sender* S;
    Receiver* R;
    
    OT(int n = N, int l = L, Sender* S = NULL, Receiver* R = NULL){
        this->n = n;
        this->l = l;
        
        GenPrime(q, 2*KAPPA);
        ZZ temp;
        do{
            RandomLen(temp, SEC_LEVEL-2*KAPPA);
            mul(temp, temp, 2);
            mul(temp, temp, q);
            add(temp, temp, 1);
        }while(!ProbPrime(temp));
        p = temp;
        sub(temp, temp, 1);
        
        find_generator(g, p, q);
        
        this->S = S;
        this->R = R;
    }
    
    void set_participants(Sender* S = NULL, Receiver* R = NULL){
        this->S = S;
        this->R = R;
        return;
    }
    
    void set_dimensions(int n, int l){
        this->n = n;
        this->l = l;
        return;
    }
    
    void first_round(){
        for(int i = 0; i < this->n; i++){
            RandomBnd(this->R->a[i], this->q);
            if(R->selection_bits[i]){
                this->S->h.x_0[i] = RandomBnd(this->p);
                this->S->h.x_1[i] = PowerMod(this->g, this->R->a[i], this->p);
            }
            else{
                this->S->h.x_0[i] = PowerMod(this->g, this->R->a[i], this->p);
                this->S->h.x_1[i] = RandomBnd(this->p);
            }
        }
        return;
    }
    
    void second_round(){
        ZZ r;
        RandomBnd(r, this->q);
        PowerMod(this->R->u, this->g, r, this->p);
        OT_Pairs keys;
        keys.x_0.SetLength(this->n);
        keys.x_1.SetLength(this->n);
        for(int i = 0; i < n; i++){
            keys.x_0[i] = PowerMod(this->S->h.x_0[i], r, this->p);
            keys.x_1[i] = PowerMod(this->S->h.x_1[i], r, this->p);
        }
        for(int i = 0; i < n; i++){
            this->R->v.x_0[i] = this->S->x.x_0[i] ^ KDF(keys.x_0[i], this->l);
            this->R->v.x_1[i] = this->S->x.x_1[i] ^ KDF(keys.x_1[i], this->l);
        }
        return;
    }
    
    void output_computation(){
        Vec<ZZ> keys;
        keys.SetLength(this->n);
        for(int i = 0; i < n; i++)
            keys[i] = PowerMod(this->R->u, this->R->a[i], this->p);
        for(int i = 0; i < n; i++)
            if(this->R->selection_bits[i])
                this->R->x[i] = this->R->v.x_1[i] ^ KDF(keys[i], this->l);
            else
                this->R->x[i] = this->R->v.x_0[i] ^ KDF(keys[i], this->l);
        return;
    }
};

class OTExtension : public OT{
public:
    
    OTExtension(int n = N, int l = L, Sender* S = NULL, Receiver* R = NULL) : OT(n, l, S, R){}
    
    void set_participants(Sender* S = NULL, Receiver* R = NULL){
        OT::set_participants(S, R);
    }
    
    void base_OT(){
        bool* selection_bits = new bool[KAPPA];
        for(int i = 0; i < KAPPA; i++)
            selection_bits[i] = RandomBits_long(1);
        this->S->R = new Receiver(false, selection_bits, KAPPA);
        OT_Pairs k;
        k.x_0.SetLength(KAPPA);
        k.x_1.SetLength(KAPPA);
        for(int i = 0; i < KAPPA; i++){
            k.x_0[i] = RandomBits_ZZ(KAPPA);
            k.x_1[i] = RandomBits_ZZ(KAPPA);
        }
        this->R->S = new Sender(false, k, KAPPA, KAPPA);
        
        int n = this->n, l = this->l;
        Sender* S = this->S;
        Receiver* R = this->R;
        
        //start = GetTime();
        OT::set_dimensions(KAPPA, KAPPA);
        OT::set_participants(this->R->S, this->S->R);
        OT::first_round();
        OT::second_round();
        OT::output_computation();
        OT::set_dimensions(n, l);
        OT::set_participants(S, R);
        //stop = GetTime();
        //cout << "BaseOT: " << (stop - start) << endl;
        return;
    }
    
    void extension_phase(){
        int offset = 0;
        int num_of_iter = (this->n/KAPPA) + 1;
        int n;
        for(int k = 0; k < num_of_iter; k++){
            if(k == num_of_iter - 1)
                n = this->n % KAPPA;
            else
                n = KAPPA;
            if(!n)
                break;
            for(int i = 0; i < KAPPA; i++)
                this->R->t[i] = PRG(R->S->x.x_0[i], n);
            // needs to be fixed to get the next KAPPA bits each time
            ZZ r;
            r = ZZ_from_bits(this->R->selection_bits, offset, n);
            start = GetTime();
            for(int i = 0; i < KAPPA; i++){
                this->S->u[i] = this->R->t[i] ^ PRG(R->S->x.x_1[i], n) ^ r;
                // needs to be fixed to get next KAPPA bits each time
                if(this->S->R->selection_bits[i])
                    this->S->q[i] = this->S->u[i] ^ PRG(this->S->R->x[i], n);
                else
                    this->S->q[i] = PRG(this->S->R->x[i], n);
            }
            stop = GetTime();
            //cout << "Compute u and q: " << (stop - start) << endl;
            start = GetTime();
            this->R->t = transpose_vector(this->R->t, n, KAPPA);
            stop = GetTime();
            //cout << "Transpose t: " << (stop - start) << endl;
            start = GetTime();
            this->S->q = transpose_vector(this->S->q, n, KAPPA);
            stop = GetTime();
            //cout << "Transpose q: " << (stop - start) << endl;
            
            ZZ s;
            s = ZZ_from_bits(this->S->R->selection_bits, 0, KAPPA);
            start = GetTime();
            for(int i = 0; i < n; i++){
                this->R->v.x_0[i] = this->S->x.x_0[i + offset] ^ CRF(this->S->q[i], this->l);
                this->R->v.x_1[i] = this->S->x.x_1[i + offset] ^ CRF(this->S->q[i] ^ s, this->l);
            }
            stop = GetTime();
            //cout << "Compute (y_0, y_1): " << (stop - start) << endl;
            
            start = GetTime();
            for(int i = 0; i < n; i++)
                if(R->selection_bits[i + offset])
                    this->R->x[i + offset] = this->R->v.x_1[i] ^ CRF(this->R->t[i], this->l);
                else
                    this->R->x[i + offset] = this->R->v.x_0[i] ^ CRF(this->R->t[i], this->l);
            stop = GetTime();
            //cout << "Compute output: " << (stop - start) << endl;
            
            offset += KAPPA;
        }
        return;
    }
};

int main(int argc, const char * argv[]) {
    double total_start, total_stop;
    start = GetTime();
    OTExtension ot;
    stop = GetTime();
    //cout << "Class Initialization: " << (stop - start) << endl;
    OT_Pairs x;
    x.x_0.SetLength(N);
    x.x_1.SetLength(N);
    bool* selection_bits = new bool[N];
    for(int i = 0; i < N; i++){
        x.x_0[i] = RandomLen_ZZ(L);
        x.x_1[i] = RandomLen_ZZ(L);
        selection_bits[i] = true*(i % 3);
    }
    Sender *S = new Sender(true, x, N, L);
    Receiver* R = new Receiver(true, selection_bits, N);
    ot.set_participants(S, R);
    start = GetTime();
    ot.base_OT();
    stop = GetTime();
    cout << "BaseOT: " << (stop - start) << endl;
    total_start = GetTime();
    ot.extension_phase();
    total_stop = GetTime();
    cout << "Extension Phase: " << (total_stop - total_start) << endl;
    /*
     ot.set_participants(S, R);
     ot.first_round();
     ot.second_round();
     ot.output_computation();
     */
    /*
     selection_bits = S->R->get_selection_bits();
     Vec<ZZ> selected_data = S->R->get_selected_data();
     x = R->S->get_data();
     */
    Vec<ZZ> selected_data = R->get_selected_data();
    ///*
    bool check_result = true;
    for(int i = 0; i < N; i++){
        if(selection_bits[i]){
            if(x.x_1[i] != selected_data[i]){
                check_result = false;
                //cout << i << " Bummer!" << endl;
            }
            else
                ;//cout << i << " Nice!" << endl;
        }
        else{
            if(x.x_0[i] != selected_data[i]){
                check_result = false;
                //cout << i << " Bummer!" << endl;
            }
            else
                ;//cout << i << " Nice!" << endl;
        }
        //cout << "Data: 0: " << x.x_0[i] << " 1: " << x.x_1[i] << endl;
        //cout << "Selection Bit: " << selection_bits[i] << " Selected Data: " << selected_data[i] << endl;
    }
    if(check_result)
        cout << "Successful Execution." << endl;
    else
        cout << "Failed Execution." << endl;
    //*/
    return 0;
}

