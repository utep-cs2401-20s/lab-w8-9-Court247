import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;

import static java.util.Arrays.*;

class AminoAcidLL{
  char aminoAcid;
  String[] codons;
  int[] counts;
  AminoAcidLL next;

  AminoAcidLL(){

  }


  /********************************************************************************************/
  /* Creates a new node, with a given amino acid/codon
   * pair and increments the codon counter for that codon.
   * NOTE: Does not check for repeats!! */

  AminoAcidLL(String inCodon){
   this.aminoAcid = AminoAcidResources.getAminoAcidFromCodon(inCodon);
   this.codons = AminoAcidResources.getCodonListForAminoAcid(aminoAcid);
   this.counts = new int [codons.length];
   this.next = null;
  } //completed

    private void addToCodons(String a){ //increases node count if codon is present
        for(int i = 0; i< counts.length; i++){
            if(codons[i]. equals(a)){
                counts[i]++;
            }
        }
    }


  /********************************************************************************************/
  /* Recursive method that increments the count for a specific codon:
   * If it should be at this node, increments it and stops,
   * if not passes the task to the next node.
   * If there is no next node, add a new node to the list that would contain the codon.
   */
  private void addCodon(String inCodon){
    if(aminoAcid == AminoAcidResources.getAminoAcidFromCodon(inCodon)){
      addToCodons(inCodon); //used helper method
    }else {
      if (next != null){
          next.addCodon(inCodon);
      }
      else{
          next = new AminoAcidLL(inCodon);
      }
    }

  } //completed


  /********************************************************************************************/
  /* Shortcut to find the total number of instances of this amino acid */
  private int totalCount(){
      int count = 0;
      for(int i =0; i < counts.length;i++){
          count += counts[i];
          }
      return count;
  }//completed

  /********************************************************************************************/
  /* helper method for finding the list difference on two matching nodes
  *  must be matching, but this is not tracked */
  private int totalDiff(AminoAcidLL inList){
    return Math.abs(totalCount() - inList.totalCount());
  }//completed


  /********************************************************************************************/
  /* helper method for finding the list difference on two matching nodes
  *  must be matching, but this is not tracked */
  private int codonDiff(AminoAcidLL inList){
    int diff = 0;
    for(int i=0; i<codons.length; i++){
      diff += Math.abs(counts[i] - inList.counts[i]);
    }
    return diff;
  }//completed

  /********************************************************************************************/
  /* Recursive method that finds the differences in **Amino Acid** counts.
   * the list *must* be sorted to use this method */
  public int aminoAcidCompare(AminoAcidLL inList){
      if(inList.isSorted() == false) {//base case
          sort(inList);
      }else{
          if (inList == null) { //if intake list is empty just return the total count of this
              int a = this.totalCount();
              return a;
          }
          if (this == null) { //if this is empty return the count for inList
              int a = this.aminoAcidCompare(inList.next);
              return a;
          }
          if (this.aminoAcid == inList.aminoAcid) {//if aminoacid is there
              return this.totalDiff(inList) + this.next.aminoAcidCompare(inList.next);
          }
          if (this.aminoAcid > inList.aminoAcid) {
              return this.totalCount() + this.next.aminoAcidCompare(inList);
          }
      }

      return inList.totalCount() + this.aminoAcidCompare(inList.next);
  }//completed

  /********************************************************************************************/
  /* Same ad above, but counts the codon usage differences
   * Must be sorted. */
  public int codonCompare(AminoAcidLL inList) { //use codondiff instead
      if(inList.isSorted() == false){ //same process for codon as it was for amino acid compare
          sort(inList);
      }else {
          if (inList == null) {
              int a = this.totalCount();
              return a;
          }
          if (this == null) {
              int a = this.aminoAcidCompare(inList.next);
              return a;
          }
          if (this.aminoAcid == inList.aminoAcid) {
              return this.codonDiff(inList) + this.next.codonCompare(inList.next);
          }
          if (this.aminoAcid > inList.aminoAcid) {
              return this.totalCount() + this.next.codonCompare(inList);
          }
      }

      return inList.totalCount() + this.codonCompare(inList.next);
  }//completed


  /********************************************************************************************/
  /* Recursively returns the total list of amino acids in the order that they are in in the linked list. */
  public char[] aminoAcidList(){
      if( next == null) {
          return new char[aminoAcid];
      }
      char[] a = next.aminoAcidList();
      char[] acidA = new char [a.length+1];
      acidA[0] = aminoAcid;
      //populate array
      for(int i = 1; i < acidA.length; i++){
          a[i] = acidA[i];
      }
      return acidA;
  }//completed

  /********************************************************************************************/
  /* Recursively returns the total counts of amino acids in the order that they are in in the linked list. */
  public int[] aminoAcidCounts(){
      if( next == null) {
          return new int[totalCount()];
      }
      int[] a = next.aminoAcidCounts();
      int[] acidA = new int [a.length+1];
      acidA[0] = totalCount();
      //populate array
      for(int i = 1; i < acidA.length; i++){
          a[i] = acidA[i];
      }
      return acidA;
  }//completed


  /********************************************************************************************/
  /* recursively determines if a linked list is sorted or not */
  public boolean isSorted(){
      AminoAcidLL head = null;
      AminoAcidLL checkSort = new AminoAcidLL();
      head = checkSort;
      if((checkSort == null) || (checkSort.next == null)) {
          return true;
      }
      return ((checkSort.aminoAcid < checkSort.next.aminoAcid) && (isSorted()));
  }//completed


  /********************************************************************************************/
  /* Static method for generating a linked list from an RNA sequence */
  public static AminoAcidLL createFromRNASequence(String inSequence){
      AminoAcidLL createNode = new AminoAcidLL();

      String s = "";
      int i = 0;
      //only get 3 characters of the string.
      for (i = 0; i < inSequence.length(); i+=3) {
          s = inSequence.substring(i, i + 3);
          //System.out.print(s + " ");
          createNode.addCodon(s.substring(0,3));
      }
      return createNode;
  }//completed


  /********************************************************************************************/
  /* sorts a list by amino acid character*/
  public static AminoAcidLL sort(AminoAcidLL inList){ //bubble sort
      AminoAcidLL sorAmino = new AminoAcidLL();
      AminoAcidLL head = sorAmino;
      AminoAcidLL index = null;
      char swap;

      if(head == null){
          return null;
      }else{
          while(head != null){
              index = head.next;
              while(index != null){
                  if(head.aminoAcid > index.aminoAcid){
                      swap = head.aminoAcid;
                      head.aminoAcid = index.aminoAcid;
                      index.aminoAcid = swap;
                  }
                  index = index.next;
              }
              head = head.next;
          }
      }
    return inList;
  }//completed

}//completed
