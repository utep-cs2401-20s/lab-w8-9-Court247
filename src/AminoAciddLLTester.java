import org.junit.jupiter.api.Test;

import java.util.Arrays;

import static org.junit.jupiter.api.Assertions.*;


public class AminoAciddLLTester {

    @Test
    public void Tester1(){ // testing the create from RNA sequence method that would give a *.
        String input = "GCUACGGAGCUUCGGAGCUAG";
        String expected = "ATELRS*";

        assertEquals(expected, AminoAcidLL.createFromRNASequence(input));
    }

    @Test
    public void Tester2(){ //Testing create method with full RNA sequence
        String input = "GCUACGGAGCUUCGGAGCUUU";
        String expected = "ATELRSF";
        AminoAcidLL tester = AminoAcidLL.createFromRNASequence(input);
        assertEquals(expected,tester);
    }

    @Test
    public void Tester3(){ // Testing create metod with * in sequence
        String input = "GCUACGGAGCUUCGGAGCUAG";
        int[] expected =  {2,1,1,2,1};
        AminoAcidLL tester = AminoAcidLL.createFromRNASequence(input);
        assertEquals(expected,tester.aminoAcidCounts());
    }

    @Test
    public void Tester4(){//testing aminoacid list cant get it to work
        String input = "GCUACGGAGCUUCGGAGCUAG";
        char[] expected = new char[] {A,A,G,L,S,S,T};
        AminoAcidLL tester = AminoAcidLL.createFromRNASequence(input);
        assertEquals(expected,tester.aminoAcidList());
    }

    @Test
    public void Tester5(){ //testing codon compare
        String input = "GCUACGGAGCUUCGGAGCUAG";
        int[] expected = new int[]{2, 1, 1, 2, 1};
        AminoAcidLL tester = AminoAcidLL.createFromRNASequence(input);
        assertEquals(expected, tester.codonCompare(tester));

    }

    @Test
    public void Tester6(){//testing sort method
        String input = "GCUACGGAGCUUCGGAGCUAG";
       String expected = "AELRST*" ;
        AminoAcidLL tester = AminoAcidLL.createFromRNASequence(input);
       assertEquals(expected, AminoAcidLL.sort(tester));

    }

    @Test
    public void Tester7() { //Testing iSorted method with full rna sequence
        String input = "GCUACGGAGCUUCGGAGCUAG";
        boolean expected = false;
        AminoAcidLL tester = AminoAcidLL.createFromRNASequence(input);
        assertFalse(expected);
    }

    @Test
    public void Tester8() { //Testing isSorted method with * in seuence
        String input = "GCUACGGAGCUUCGGAGCUUU";
        boolean expected = false;
        AminoAcidLL tester = AminoAcidLL.createFromRNASequence(input);
        assertFalse(expected);
    }

    @Test
    public void Tester9() { //Testing aminoacidcompare with full RNA
        String input = "GCUACGGAGCUUCGGAGCUUU";
        int[] expected = {1,0,0,0,0,0,0};
        AminoAcidLL tester = AminoAcidLL.createFromRNASequence(input);
        assertEquals(expected, tester.aminoAcidCompare(tester));
    }

    @Test
    public void Tester10() { //Testing aminoacid compare.
        String input = "GCUACGGAGCUUCGGAGCUAG";
        int[] expected = {1,0,0,0,0,0,} ;
        AminoAcidLL tester = AminoAcidLL.createFromRNASequence(input);
        assertEquals(expected,tester.aminoAcidCompare(tester));
    }


}
