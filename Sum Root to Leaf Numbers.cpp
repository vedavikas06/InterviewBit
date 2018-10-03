/**
 * Definition for binary tree
 * struct TreeNode {
 *     int val;
 *     TreeNode *left;
 *     TreeNode *right;
 *     TreeNode(int x) : val(x), left(NULL), right(NULL) {}
 * };
 */
 
 int sum(TreeNode *A,int val){
     
     if(A==NULL){return 0;}
     if(!A->left && !A->right){
         return (val*10+A->val)%1003;
     }
     
     return (sum(A->left,(val*10+A->val)%1003)+sum(A->right,(val*10+A->val)%1003))%1003;
     
 }
int Solution::sumNumbers(TreeNode* A) {
    
    return sum(A,0);
}
