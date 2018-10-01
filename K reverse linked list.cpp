/**
 * Definition for singly-linked list.
 * struct ListNode {
 *     int val;
 *     ListNode *next;
 *     ListNode(int x) : val(x), next(NULL) {}
 * };
 */
ListNode* Solution::reverseList(ListNode* A, int B) {
    
        ListNode *current = A;
        ListNode *prev = NULL, *next = NULL,*last=NULL,*last1=NULL;
        int i =0,j=0;
 
        if(B!=1){
        while (i<B)
        {
            if(i==0){
                last1=current;
            }
            
            if(last){
                last->next = current;
            }
            // Store next
            next = current->next;
 
            // Reverse current node's pointer
            current->next = prev;
 
            // Move pointers one position ahead.
            prev = current;
            current = next;
            
            i++;
            
            if(i==B){
                if(j==0){
                    A = prev;
                    j++;
                }
                prev = NULL;
                i=0;
                last = last1;
            }
        }
        
        }
        return A;
}
