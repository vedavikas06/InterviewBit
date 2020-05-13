class Solution {
public:
    string removeKdigits(string num, int k) {
        int len = num.length();
        if(k==len){
            return "0";
        }
        string small = "",pick;
        int no_times = len-k;
        for(int i=0;i<no_times;i++){
            pick = num.substr(0,k+1);
            int minn = 0;
            for(int j=1;j<pick.length();j++){
                if(pick[j]<pick[minn]){
                    minn = j;
                }
            }
            small += pick[minn];
            if(minn+1<num.length()){
                num = num.substr(minn+1);    
            }
            
            k=k-minn;
            
        }
        
        
        int last = -1;
        for(int i=0;i<small.length();i++){
            if(small[i]=='0' && last == i-1){
                last = i;
            }else{
                break;
            }
        }
        if(last==small.length()-1){
            return "0";
        }
        return small.substr(last+1);
        
    }
};
