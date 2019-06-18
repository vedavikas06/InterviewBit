//BFS kind of problem
vector<int> Solution::solve(int A, int B, int C, int D) {
    
    priority_queue<int,vector<int> ,greater<int> > p;
    
    vector<int> v;
    v.push_back(A);v.push_back(B);v.push_back(C);
    for(int i =0;i<v.size();i++)
    p.push(v[i]);
    
    int cnt=0;
    map<int,int> m;
    vector<int> res;
    while(cnt < D){
        int tp = p.top();
        p.pop();
        
        if(m.find(tp)==m.end()){
        res.push_back(tp);
        
        for(int i =0;i<v.size();i++){
            p.push(v[i]*tp);
        }
        
        cnt++;
        m[tp] = 1;
        }
        
        
        
    }
     return res; 
    
}
