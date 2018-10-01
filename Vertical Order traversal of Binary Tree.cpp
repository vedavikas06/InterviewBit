/**
 * Definition for binary tree
 * struct TreeNode {
 *     int val;
 *     TreeNode *left;
 *     TreeNode *right;
 *     TreeNode(int x) : val(x), left(NULL), right(NULL) {}
 * };
 */
vector<vector<int> > Solution::verticalOrderTraversal(TreeNode* A) {
    
    map<int,vector<int>> m;
    
    queue<pair<TreeNode *,int> > qu;
    
     vector<vector<int> > vec;
     
     if(A){
    
    qu.push({A,0});
    
    m[0].push_back(A->val);
    
   
    
    
    while(!qu.empty()){
       pair<TreeNode *,int> ext = qu.front();
       qu.pop();
       
       if(ext.first->left){
           qu.push({ext.first->left,ext.second-1});
           m[ext.second-1].push_back(ext.first->left->val);
       }
       
       if(ext.first->right){
           qu.push({ext.first->right,ext.second+1});
           m[ext.second+1].push_back(ext.first->right->val);
       }
        
        
        
        
    }
    
    for(auto i:m){
       vec.push_back(i.second); 
    }
    
    return vec;
    
     }
     
     return vec;
    
}
