class Solution {
public:
    bool static temp(vector<int> a,vector<int> b){
        return a[0] < b[0];
    }
    vector<vector<int>> merge(vector<vector<int>>& intervals) {
        int n = intervals.size();
        sort(intervals.begin(),intervals.end(),temp);
        
        stack<vector<int>> s;
        
        for(int i=0;i<n;i++){
            if(s.empty()){
                s.push(intervals[i]);
            }else{
                vector<int> pr = intervals[i];
                while(!s.empty() && s.top()[1] >= pr[0]){
                    vector<int> tp = s.top();
                    s.pop();
                    pr[0] = tp[0];
                    pr[1] = max(tp[1],pr[1]);
                }
                s.push(pr);
            }
        }
        vector<vector<int>> op;
        while(!s.empty()){
            op.insert(op.begin(),s.top());
            s.pop();
            
        }
        return op;
    }
};
