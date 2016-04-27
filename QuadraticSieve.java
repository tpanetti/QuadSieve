import java.util.Scanner;
import java.util.ArrayList;
import java.util.HashSet;
/*
**This program will run a quadratic sieve to
**factor some numbers
**@author Thomas Panetti & Noemi Glaeser
**
**
*/
 
public class QuadraticSieve
{
  //This is a private class
  //to be used as a touple for
  //saving an ArrayList and it's integer index
  private static final int REDUCERANGE = 8;
  private static class Pair<X, Y> 
  {
    private final X x;
    private final Y y;
    public Pair(X x, Y y) 
    {
      this.x = x;
      this.y = y;
    }
    public X getX(){ return x;}
    public Y getY(){ return y;}
  }


  public static void main(String[] args)
  {
    //Grab number N to be factored
    long n;
    long factorbase;
    //take from command line
    if(args.length == 2)
    {
      n = Long.parseLong(args[0]);
      factorbase = Long.parseLong(args[1]);
    } 
    //else take from scanner
    else
    {
      Scanner input = new Scanner(System.in);
      System.out.println("Enter an N: ");
      n = input.nextLong();
	  System.out.println("Enter a factorbase: ");
      factorbase = input.nextLong();
    }
    //Find R
    long R = (long)Math.sqrt(n);

    //generate a list of primes up to factorbase
    ArrayList<Integer> primes = eratosthenesSieve(factorbase);
           //for(int i = 0; i < primes.size(); i++)
            //System.out.println(primes.get(i) + " "); 
    //generate a list that satisfies the eqn 
    ArrayList<Long> bSmooth = findSmoothness(R/REDUCERANGE, n);
    ArrayList<Long> residues = calcResiduals(primes, bSmooth);
  
    //Refactor and Gauss
    ArrayList<Pair<ArrayList<Long>, Integer>> refactoredResidual = refactor(residues, bSmooth, primes);
    //reduce the matrix mod 2
    ArrayList<Pair<ArrayList<Integer>, Integer>> reducedResidual = reduceModTwo(refactoredResidual);
    
    //print reduced residual matrix
    //System.out.println("Reduced matrix is: ");
    //int i = 0;
    //for(Pair<ArrayList<Integer>,Integer> pair : reducedResidual)
    //{
      //System.out.print(i + " ");
      //System.out.print("[");
      //for(long j : pair.getX())
        //System.out.print(j + ", ");
      //System.out.println("]");
    //i++;
    //}
    //perform Gauss-Jordan elimination on the matrix
    //build list of hits I.E. (Index, r1, r2, r3, etc)
    ArrayList<ArrayList<Integer>> gaussHits = gauss(reducedResidual);
    //rebuild equation
    HashSet<Long> factors = rebuildEquation(gaussHits, n);
    System.out.print("The factors are: " );   
    for(long factor : factors)
      System.out.print(factor + ", ");
    
  }
  
  /**
  * This method will perform an eratostheness Sieve to calculate all primes
  * under the input factobase
  *
  * @param size This is the size to create the array of primes
  * @return 	An ArrayList of primes under size 
  */        
  public static ArrayList<Integer> eratosthenesSieve(long size) 
  {
    ArrayList<Integer> primes = new ArrayList<Integer>();
    //populate a list with i = 2 to n
    for(int i = 2; i < size; i++) 
      primes.add(i);
    
    //set p equal to 2?
    int p = 2;
    //loop
    while(Math.pow(p, 2) < size)
    {
      for(int i = primes.indexOf(p) + 1; i < primes.size(); i++)
      {
        ////System.out.println(primes.get(i) +": divided by " + p);
        if(primes.get(i) % p == 0)
        {
          primes.remove(i);
          
        }
      }
    p++;
    }
    primes.add(0, -1);
    return primes;
  }
  
  /**
  *
  *@param numToFactor 	the input number we are factoring
  *@param range	the range at which to create the matrix from -range to range
  *@return 	An ArrayList of longs representing the smooth values found
  */
  public static ArrayList<Long> findSmoothness(long range, long numToFactor)
  {
    //Offset = R. Change later
    ArrayList<Long> Qs = new ArrayList<Long>();
    //System.out.println("Root is " + Math.sqrt(numToFactor) + "\tindex\t" + "(r+n)^2\t" + "(r+n)^2 - N)");
    for(long i = -range; i <= range; i++) 
    {
      //store log of Q possibly
      long Q = ((long)(Math.pow((Math.sqrt(numToFactor) + i), 2)));
       
      long save = Q - numToFactor;
      Qs.add(save);
      //System.out.println(i +"\t" + (Qs.size() - 1) + "\t" + Q + "\t" + save);
    }
  return Qs;
  }
  
  public static ArrayList<Long> calcResiduals(ArrayList<Integer> primes, ArrayList<Long> smoothlist)
  {
    //ArrayList<Pair<Integer, Integer>> residuals = new ArrayList<Pair<Integer,Integer>>();

    ArrayList<Long> copy = new ArrayList<Long>();
    for(long i : smoothlist)
      copy.add(i);
    int start;
    for(int i : primes)
    {
    if(i == -1)
      continue;
    //System.out.println("Sieve with " + i);
      for(int j = 0; j < copy.size(); j++)
      {
        if(copy.get(j) % i == 0)
	{
          start = j;
          int index = j;
          do{ 
            long temp = copy.get(j);
            while(temp != 0 && temp % i == 0)
            {
              //System.out.println("Temp: " + temp + "\ti: " + i + "temp % i = " + (temp%i));
              //System.out.println("Divide\t" + i + "\t " + (j-(smoothlist.size()/2)) + "\t " + temp + "\t" + (temp / i));
              temp = temp / i;
            }
            copy.set(j, temp);
            j = (j + i) % copy.size();
          } while(j != start);
	}  
      }
      //System.out.print("Array [");
      /* for(long num : copy)
      {
        //System.out.print(num + ", "); 
      } */
      //System.out.println("]");
    }
  //Take abs value of the array
  for(int i = 0; i < copy.size(); i++)
    copy.set(i, Math.abs(copy.get(i)));
  return copy;
  
  }


  public static ArrayList<Pair<ArrayList<Long>,Integer>> refactor(ArrayList<Long> residues, ArrayList<Long> original,
							ArrayList<Integer> primes)
  {
    //initialize arraylist and fill with zero arryalists
    ArrayList<Pair<ArrayList<Long>, Integer>> exponents = new ArrayList<Pair<ArrayList<Long>, Integer>>();
    //cheating
    long zero = 0;
    for(int i = 0; i < residues.size(); i++)
    {
      if(residues.get(i) == 1)
      {
        ArrayList<Long> exponent = new ArrayList<Long>();
        for(int index = 0; index < primes.size(); index++)
          exponent.add(zero);
	//System.out.println("Smooth\t" + (i - (original.size()/2)) + "\t" + i + "\t" + original.get(i));
        long temp = original.get(i);
	if(temp < 0)
	{
	  temp = temp * -1;
	  exponent.set(0,new Long(1));
	}
        for(int pIndex = 1; pIndex < primes.size(); pIndex++) //pIndex set to 1 to skip the "prime" -1 
        {
          while(temp % primes.get(pIndex) == 0)
	        {
	          //System.out.println("DIVIDE\t" + pIndex +"\t" + primes.get(pIndex) + "\t" + temp);
	          temp = temp / primes.get(pIndex);
	          exponent.set(pIndex, (exponent.get(pIndex))+1);
	        }
	
	      }     
        exponents.add(new Pair<ArrayList<Long>, Integer>(exponent, i));
        }
       
      //System.out.print("[");
      /* if(exponents.size() > 0)
      {
        for(long whatever : exponents.get(exponents.size() -1).getX())
        {
          //System.out.print(whatever +", ");
        }
        //System.out.println("]");
      } */
    }


    return exponents;
  }
  public static ArrayList<Pair<ArrayList<Integer>, Integer>> reduceModTwo(ArrayList<Pair<ArrayList<Long>, Integer>> matrix)
  {
    ArrayList<Pair<ArrayList<Integer>, Integer>> modTwo = new ArrayList<Pair<ArrayList<Integer>, Integer>>();
    for(int i = 0; i < matrix.size(); i++)
    {
      //System.out.print("Row: " + i + "\t[");  
      
      for(int j = 0; j < matrix.get(i).getX().size(); j++)
      {
        //System.out.print(matrix.get(i).getX().get(j) + ", ");
        matrix.get(i).getX().set(j, matrix.get(i).getX().get(j) % 2);
  
      }
      //System.out.print("]\nReduced\t[");
      //for(long temp : matrix.get(i).getX())
        //System.out.print(temp + ", ");
      //System.out.println("]");
     }
     for(Pair<ArrayList<Long>, Integer> pair : matrix)
     {
       ArrayList<Integer> ints = new ArrayList<Integer>();
       for(long number : pair.getX())
         ints.add((int)number); 
       modTwo.add(new Pair<ArrayList<Integer>, Integer>(ints, pair.getY()));
     }
     return modTwo;
  }

  
  
  public static ArrayList<ArrayList<Integer>> gauss(ArrayList<Pair<ArrayList<Integer>,Integer>> array)
  {
    ArrayList<Pair<ArrayList<Integer>, ArrayList<Integer>>> reducedMatrix = new ArrayList<Pair<ArrayList<Integer>, ArrayList<Integer>>>();
    //initialize ArrayList
    for(int i = 0; i < array.size(); i++)
    {
      ArrayList<Integer> rowReductions = new ArrayList<Integer>();
      //first value in the index is the original index of said row
      rowReductions.add(array.get(i).getY());
      reducedMatrix.add(new Pair<ArrayList<Integer>, ArrayList<Integer>>(array.get(i).getX(), rowReductions));
    }
    int j = 0;
    for(int i = 0; i < reducedMatrix.size(); i++)
    { 
      if(j >= reducedMatrix.get(i).getX().size())
        break;
      if(reducedMatrix.get(i).getX().get(j) == 0)
      {
        //index of the next row that has a 1 to swap for the zero
        int oneRow;
        for(oneRow = i + 1; oneRow < reducedMatrix.size(); oneRow++)
        {
          //break at 1 to preserve oneRow for swapping
          if(reducedMatrix.get(oneRow).getX().get(j) == 1)
            break;
        }
        //if it's not equal to the size, we did find a 1
        if(oneRow < reducedMatrix.size())
        {  //one found, swap rows
          Pair<ArrayList<Integer>, ArrayList<Integer>> temp = reducedMatrix.get(i);
          reducedMatrix.set(i, reducedMatrix.get(oneRow));
          reducedMatrix.set(oneRow, temp);
        }
        //else move to the next column
        else
        {
          i--;
          j++;
          continue;
        }
      }
      
      //else we have 1, zero out the column
      for(int row = i + 1; row < reducedMatrix.size(); row++)
      {
        if(reducedMatrix.get(row).getX().get(j) == 0)
          continue;
        reducedMatrix.get(row).getX().set(j,0);
        reducedMatrix.get(row).getY().add(i);
      }
      
      j++;
    }

    ArrayList<ArrayList<Integer>> zeroSums = new ArrayList<ArrayList<Integer>>(); 

    for(int i = 0; i < reducedMatrix.size(); i++)
    {
      for(j = 0; j < reducedMatrix.get(i).getX().size(); j++)
      {
        if(reducedMatrix.get(i).getX().get(j) != 0)
        { 
          //break out of the loop by ending condition
          j = reducedMatrix.get(i).getX().size();
          continue;
        }
      }
      zeroSums.add(reducedMatrix.get(i).getY());
    }  
    return zeroSums;   
  }
  
  /**
  * This method will rebuild the equation x^2 = y^2
  *
  *
  *
  */
  public static HashSet<Long> rebuildEquation(ArrayList<ArrayList<Integer>> gaussHits, long bigN)
  {
    int R = ((int)(Math.sqrt(bigN))) / REDUCERANGE;
    HashSet<Long> factors = new HashSet<Long>();
    for(ArrayList<Integer> hit : gaussHits)
    {
      int lhs = 1;
      int rhs = 1;
      for(int num = 1; num < hit.size(); num++)
      {
        lhs *= Math.pow(R + (hit.get(num)- R), 2);
        rhs *= Math.pow(R + (hit.get(num) - R), 2) - bigN;
      }
      long x = (long)(Math.sqrt(Math.abs(lhs)));
      long y = (long)(Math.sqrt(Math.abs(rhs)));
      if(x == y)
        continue;
      long factor1 = (gcd(bigN, (x+y)));
      long factor2 = (gcd(bigN, (x-y)));
      ////System.out.print("Hit\t[");
      //for(long number : hit)
      //  //System.out.print(number + ", "); 
      //System.out.println("]\n LHS: " + lhs + ".\t RHS: " + rhs);
      //System.out.println("X: " + x + "\tY: " + y);
      //System.out.println("gcd(" + bigN +", " + (x+y) +": " + factor1 +
                        //"\ngcd(" + bigN + ", " + (x-y) + ": " + factor2);
      
      factors.add(Math.abs(factor1));
      factors.add(Math.abs(factor2));
      //check if the factors make other factors!
      if(bigN % factor1 == 0)
        factors.add(Math.abs(bigN/factor1));
      if(bigN % factor2 == 0)
        factors.add(Math.abs(bigN/factor2));
       if(factors.size() > 1)
         factors = completeFactors(factors);
    }
      
    return factors;
  }
  
  public static HashSet<Long> completeFactors(HashSet<Long> factors) 
  {
    HashSet<Long> newFactors = new HashSet<Long>();
    for(long fact : factors)
      newFactors.add(fact);
    for(long fact : factors)
      for(long fact2 : factors)
        if(fact > fact2 && fact % fact2 == 0)
          newFactors.add(fact / fact2);
    return newFactors;
  }
  
  /**
  * This method will run the euclidean algorithm to
  * factor two integer into their greatest common denominator
  *
  *@param a 	An integer to be reduced
  *@param b 	An integer to be reduced
  *@return	The greatest common denominator between a and b
  */
  public static long gcd(long a, long b)
  {
    while(b != 0)
    {
      long t = b;
      b = a % b;
      a = t;
    }
    return a;
  }



}
