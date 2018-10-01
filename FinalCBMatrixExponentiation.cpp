//Codersbit Final question on Matrix Exponentiation
#include <bits/stdc++.h>

using namespace std;

#define pb push_back
#define sf scanf
#define pf printf
#define f first
//#define s second
#define clr(x,y) memset(x,y,sizeof x)
#define LL long long

#define mx 100009


long long mod = 1000000007;
// #define LL long long
vector<vector<LL> > multiply(vector<vector<LL> >g1, vector<vector<LL> >g2)
{

    vector<vector<LL> >mul(g1.size(), vector<LL> (g1.size(), 0));
    for (LL i = 1; i < g1.size(); i++)
    {
        for (LL j = 1; j < g1.size(); j++)
        {
            mul[i][j] = 0;
            for (LL k = 1; k < g1.size(); k++)
                mul[i][j] = (mul[i][j] + g1[i][k] * g2[k][j]) % mod;
        }
    }
    return mul;

}


LL solve(LL A, LL B) {

    vector<LL> cnt(B + 1, 0);
    vector<vector<LL> > tot(B + 1, vector<LL>(B + 1, 0));
    vector<vector<LL> > fin(B + 1, vector<LL>(B + 1, 0));
    for (LL i = 1; i <= B; i++) {
        for (LL j = 1; j <= B; j++) {
            if (__gcd(i, j) == 1) {
                cnt[i]++;
                tot[i][j] = 1;
            }
            
           // cout << tot[i][j] << " ";
        }
        fin[i][i] =1;

       // cout << endl;
    }


    vector<LL> bin;

    LL n = A-1;

    // A==1 i.e ans is B then
    if(n==0){
     return B;   
    }

    while (n != 0) {
        bin.pb(n % 2);
        n = n / 2;
    }

    // for(auto i:bin)
    //     cout << i << " ";
    // cout << endl;
    
    vector<vector<LL> > g = tot;
    

    for (LL i = 0; i < bin.size(); i++) {
        if (bin[i] != 0 ) {
            
            fin = multiply(fin,g);
    

        }
        g = multiply(g, g);

    }

   

    long long val = 0;
     for (LL i = 0; i < g.size(); i++)
                {
                    for (LL j = 0; j < g.size(); j++)
                    {
                       
                            val = (val + fin[i][j] ) % mod;
                    }
                }
        

          


    return val;

}


  main() {

    // LL start_s = clock();

    ios_base::sync_with_stdio(false);

    cin.tie(NULL);

    LL a, b;

    cin >> a >> b;

    cout << solve(a, b) << endl;




    // LL stop_s = clock();

    // cout << "time: " << (stop_s-start_s)/double(CLOCKS_PER_SEC) << endl;

    return 0;


}


