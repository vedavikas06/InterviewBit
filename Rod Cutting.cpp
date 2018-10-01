void simplify(int l,int r,map<pair<int,int>,long long> &dp,vector<int> &B,map<pair<int,int>,int> &st){
    
    
    
    if(r-l <= 1){
        dp[{l,r}] = 0;
    }
    
    if(dp.find({l,r})!=dp.end()){
        return ;
    }
    // if(r-l ==2){
    //     vector<int> vh;
    //     vh.push_back(B[l+1]);
    //     return {B[r]-B[l],vh};
    // }
    long long minn = INT_MAX,mn=l+1;
    
    for(int i=l+1;i<r;i++){
        // pair<int,vector<int>> a = simplify(l,i,dp,B),b=simplify(i,r,dp,B);
        if(dp.find({l,i})==dp.end()){
        simplify(l,i,dp,B,st);
        }
        
        if(dp.find({i,r})==dp.end()){
        simplify(i,r,dp,B,st);
        }
        
    
    
        int val = B[r]-B[l] + dp[{l,i}] + dp[{i,r}];
        
        if(val < minn){
        
        
        minn = val;
         
        mn = i;
         
        }
    
        
    }
    
     dp[{l,r}] = minn;
         
     st[{l,r}] = mn;
    
    
    
    return ;
    
}


void extract(int l,int r,map<pair<int,int>,int> st,vector<int> &B,vector<int> &ans){
    if(r-l<=1){
        return;
    }
    
    int inx = st[{l,r}];
    
    ans.push_back(B[inx]);
    extract(l,inx,st,B,ans);
    extract(inx,r,st,B,ans);
    
    
    
}


vector<int> Solution::rodCut(int A, vector<int> &B) {

    
    sort(B.begin(),B.end());
    
    B.insert(B.begin(),0);
    
    B.insert(B.end(),A);
    
    // for(auto i:B){
    //     cout << i << " ";
    // }
    // cout << endl;
    
    
    map<pair<int,int>,long long> dp;
    map<pair<int,int>,int> st;

    simplify(0,B.size()-1,dp,B,st);
    // cout << val.first << endl;
    // for(auto i:val.second){
    //     cout << i << " ";
    // }
    // cout << endl;
    
    
    vector<int> ans;
    
    extract(0,B.size()-1,st,B,ans);
    
    return ans;
    
}


// void solve(vector<int> &A,int l,int r,vector<int> &ans,map< pair<int,int> , long long > &mp,map< pair<int,int> , int > &index){
                
//     if(r-l<=1) {
//         mp[{l,r}]=0;
//     }
    
//     if(mp.find({l,r})!=mp.end()) return;
    
//     long long cost=1e15;
//     int mn=l+1;
    
//     for(int i=l+1;i<r;i++){
//         if(mp.find({l,i})==mp.end()) solve(A,l,i,ans,mp,index);
//         if(mp.find({i,r})==mp.end()) solve(A,i,r,ans,mp,index);
        
//         if(cost > mp[{l,i}] + mp[{i,r}] + A[r]-A[l]){
//             cost=mp[{l,i}] + mp[{i,r}] + A[r]-A[l];
//             mn=i;
//         }
//     }
//     mp[make_pair(l,r)]=cost;
//     index[{l,r}]=mn;
    
//     //cout<<l<<" "<<r<<" "<<mn<<" " <<cost+r-l <<endl;
// }

// void fillA(int l,int r,vector<int> &A,vector<int> &ans,map< pair<int,int> , long long > &mp,map< pair<int,int> , int > &index){
//     if(r-l<=1) return;
    
//     int i=index[{l,r}];
//     ans.push_back(A[i]);
    
//     //cout<<l<<" "<<r<<" "<<i<<endl;
    
//     fillA(l,i,A,ans,mp,index);
//     fillA(i,r,A,ans,mp,index);
// }

// vector<int> Solution::rodCut(int N, vector<int> &A) {
//     vector<int> ans,B;
    
//     B.push_back(0);
//     for(auto x:A) B.push_back(x);
//     B.push_back(N);
    
//     int n=B.size();
    
//     map< pair<int,int> , long long > mp;
//     map< pair<int,int> , int > index;
    
//     solve(B,0,n-1,ans,mp,index);
    
//     fillA(0,n-1,B,ans,mp,index);
    
//     return ans;
// }


