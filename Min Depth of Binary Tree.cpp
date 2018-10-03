/**
 * Definition for binary tree
 * struct TreeNode {
 *     int val;
 *     TreeNode *left;
 *     TreeNode *right;
 *     TreeNode(int x) : val(x), left(NULL), right(NULL) {}
 * };
 */
int Solution::minDepth(TreeNode* A) {
    if(A==NULL){return 0;}
    else if(!A->left && !A->right){return 1;}
    
    if(!A->left){
        return minDepth(A->right)+1;
    }else if (!A->right){
        return minDepth(A->left)+1;
    }
    return 1+min(minDepth(A->left),minDepth(A->right));
}
