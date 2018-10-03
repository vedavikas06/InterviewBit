/**
 * Definition for binary tree
 * struct TreeNode {
 *     int val;
 *     TreeNode *left;
 *     TreeNode *right;
 *     TreeNode(int x) : val(x), left(NULL), right(NULL) {}
 * };
 */
 
 vector<vector<int> > b;
 
 
void sum(TreeNode* A,vector<int> v, int B){
    
  if(A==NULL){
      return;
  }  
  if(!A->left && !A->right ){
      if(B==A->val){
          v.push_back(A->val);
          b.push_back(v);
          
      }
      return;
  }
   v.push_back(A->val);
   sum(A->left,v,B-(A->val)); 
   sum(A->right,v,B-(A->val));
   
}
vector<vector<int> > Solution::pathSum(TreeNode* A, int B) {
    
    vector<int> v;
    sum(A,v,B);
    
    vector<vector<int> > cpy = b;
    
    b.clear();
    
    return cpy;
}
