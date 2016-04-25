import java.util.Scanner;
import java.util.ArrayList;
import java.util.HashSet;
/*
**This program will run a quadratic sieve to find
**some numbers? Definitely
**@author Thomas Panetti & Noemi
**
**
*/
 
  //Numbers that work
  //1037, 100, 15, 45
  //Numbers that don't work
  //most, 437, 1147, 2491 
public class QuadraticSieve
{

  //private final int SMOOTHNESSBOUND = 1000;
  //DELETE offset (No public vars) 
  //private static long offset;
  //This is a private class
  //to be used as a touple for
  //saving an ArrayList and it's integer index
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
    //Grab number N
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
      factorbase = input.nextLong();
    }
    
    //Find R?
    long R = (long)Math.sqrt(n);  //as long as the sqrt is an integer

    //generate a list of "Big Roots"
    ArrayList<Integer> primes = eratosthenesSieve(factorbase);
           for(int i = 0; i < primes.size(); i++)
            System.out.println(primes.get(i) + " "); 
    //generate a list that satisfies the eqn 
    ArrayList<Long> bSmooth = findSmoothness(20, n); //MAGIC NUMBER 20
    ArrayList<Long> residues = calcResiduals(primes, bSmooth);
  
    //Refactor and Gauss
    ArrayList<Pair<ArrayList<Long>, Integer>> refactoredResidual = refactor(residues, bSmooth, primes);
    //reduce the matrix mod 2
    ArrayList<Pair<ArrayList<Integer>, Integer>> reducedResidual = reduceModTwo(refactoredResidual);
    
    //print reduced Residual matrix
    System.out.println("Reduced matrix is: ");
    int i = 0;
    for(Pair<ArrayList<Integer>,Integer> pair : reducedResidual)
    {
      System.out.print(i + " ");
      System.out.print("[");   
      for(long j : pair.getX())
        System.out.print(j + ", ");
      System.out.println("]");
    i++;
    }
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
  * This method will perform a Sieve to calculate all primes
  * under the input 
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
        //System.out.println(primes.get(i) +": divided by " + p);
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
  
  public static ArrayList<Long> findSmoothness(long range, long numToFactor)
  {
    //Offset = R. Change later
    ArrayList<Long> Qs = new ArrayList<Long>();
    System.out.println("Root is " + Math.sqrt(numToFactor) + "\tindex\t" + "(r+n)^2\t" + "(r+n)^2 - N)");
    for(long i = -range; i <= range; i++) 
    {
      //store log of Q possibly
      long Q = ((long)(Math.pow((Math.sqrt(numToFactor) + i), 2)));
       
      long save = Q - numToFactor;
      Qs.add(save);
      System.out.println(i +"\t" + (Qs.size() - 1) + "\t" + Q + "\t" + save);
    }
  return Qs;
  }
  
  public static ArrayList<Long> calcResiduals(ArrayList<Integer> primes, ArrayList<Long> smoothlist)
  {
    //ArrayList<Pair<Integer, Integer>> residuals = new ArrayList<Pair<Integer,Integer>>();

    ArrayList<Long> copy = new ArrayList<Long>();
      System.out.println("Size of smoothlist should be 41, size is: " + smoothlist.size());
    for(long i : smoothlist)
      copy.add(i);
    int start;
    for(int i : primes)
    {
    if(i == -1)
      continue;
    System.out.println("Sieve with " + i);
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
              System.out.println("Temp: " + temp + "\ti: " + i + "temp % i = " + (temp%i));
              System.out.println("Divide\t" + i + "\t " + (j-20) + "\t " + temp + "\t" + (temp / i));//MAGIC NUMBER
              temp = temp / i;
            }
            copy.set(j, temp);
            j = (j + i) % copy.size();
          } while(j != start);
	}  
      }
      System.out.print("Array [");
      for(long num : copy)
      {
        System.out.print(num + ", "); 
      }
      System.out.println("]");
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
	System.out.println("Smooth\t" + (i - 20) + "\t" + i + "\t" + original.get(i)); //MAGIC NUMBER
        long temp = original.get(i);
	if(temp < 0)
	{
	  temp = temp * -1;
	  exponent.set(0,new Long(1));
	}
        for(int pIndex = 1; pIndex < primes.size(); pIndex++) //pIndex set to one to skip -1 
        {
          while(temp % primes.get(pIndex) == 0)
	        {
	          System.out.println("DIVIDE\t" + pIndex +"\t" + primes.get(pIndex) + "\t" + temp);
	          temp = temp / primes.get(pIndex);
	          exponent.set(pIndex, (exponent.get(pIndex))+1);
	        }
	
	      }     
        exponents.add(new Pair<ArrayList<Long>, Integer>(exponent, i));
        }
       
      System.out.print("[");
      if(exponents.size() > 0)
      {
        for(long whatever : exponents.get(exponents.size() -1).getX())
        {
          System.out.print(whatever +", ");
        }
        System.out.println("]");
      }
    }


    return exponents;
  }
  public static ArrayList<Pair<ArrayList<Integer>, Integer>> reduceModTwo(ArrayList<Pair<ArrayList<Long>, Integer>> matrix)
  {
   //TODO Print out before and after reduction
   //smooth -1 19 is wrong double check  
    ArrayList<Pair<ArrayList<Integer>, Integer>> modTwo = new ArrayList<Pair<ArrayList<Integer>, Integer>>();
    for(int i = 0; i < matrix.size(); i++)
    {
      System.out.print("Row: " + i + "\t[");  
      
      for(int j = 0; j < matrix.get(i).getX().size(); j++)
      {
        System.out.print(matrix.get(i).getX().get(j) + ", ");
        matrix.get(i).getX().set(j, matrix.get(i).getX().get(j) % 2);
  
      }
      System.out.print("]\nReduced\t[");
      for(long temp : matrix.get(i).getX())
        System.out.print(temp + ", ");
      System.out.println("]");
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
        //System.out.println("i = " + i);
        //System.out.println("j = " + j);
        for(oneRow = i + 1; oneRow < reducedMatrix.size(); oneRow++)
        {
          //System.out.println("oneRow = " + oneRow);
          //System.out.println("matrix size = " + reducedMatrix.size());
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
    //Old way
    /*System.out.println("Enter Gauss"); 
    ArrayList<ArrayList<Integer>> gaussHits = new ArrayList<ArrayList<Integer>>();
    //calc doubles
    for(int i = 0; i < array.size(); i++)
    {
      for(int j = i+1; j < array.size(); j++)
      {
        for(int index = 0; index < array.get(i).getX().size(); index++)
        {
	  if((array.get(i).getX().get(index) + array.get(j).getX().get(index)) % 2 != 0)
	    break;
          //BAD CODE
	  if(index == array.get(i).getX().size() - 1) //THIS MEANS WERE AT THE END
	  {
            System.out.println("Hit\t" + (i) + "\t" + (j)); 
 	    ArrayList<Integer> thisHit = new ArrayList<Integer>();
            thisHit.add(i);
            thisHit.add(j);
            gaussHits.add(thisHit);
          }
        }
  
      } 
    }
    //calc triples
    for(int i = 0; i < array.size(); i++)
    {
      for(int j = i+1; j < array.size(); j++)
      {
        for(int k = j + 1; k < array.size(); k++)
        {
          for(int index = 0; index < array.get(i).getX().size(); index++)
          {
	    if((array.get(i).getX().get(index) + array.get(j).getX().get(index) + array.get(k).getX().get(index)) % 2 != 0)
	      break;
            //BAD CODE
	    if(index == array.get(i).getX().size() - 1) //THIS MEANS WERE AT THE END
	    {
              System.out.println("Hit\t" + (i) + "\t" + (j) + "\t" + (k)); 
              ArrayList<Integer> thisHit = new ArrayList<Integer>();
              thisHit.add(i);
              thisHit.add(j);
              thisHit.add(k);
              gaussHits.add(thisHit);	      
            }
  	  }
        }
      } 
    }
    return gaussHits;
    */
  }
  
  /**
  * This method will rebuild the equation (x
  *
  *
  *
  */
  public static HashSet<Long> rebuildEquation(ArrayList<ArrayList<Integer>> gaussHits, long bigN)
  {
    int R = (int)(Math.sqrt(bigN));
    HashSet<Long> factors = new HashSet<Long>();
    for(ArrayList<Integer> hit : gaussHits)
    {
      int lhs = 1;
      int rhs = 1;
      for(int num = 1; num < hit.size(); num++)
      {
        lhs *= Math.pow(R + (hit.get(num) - 20), 2);
        rhs *= Math.pow(R + (hit.get(num) - 20), 2) - bigN;
      }
      long x = (long)(Math.sqrt(Math.abs(lhs)));
      long y = (long)(Math.sqrt(Math.abs(rhs)));
      if(x == y)
        continue;
      long factor1 = (gcd(bigN, (x+y)));
      long factor2 = (gcd(bigN, (x-y)));
      System.out.print("Hit\t[");
      for(long number : hit)
      	System.out.print(number + ", "); 
      System.out.println("]\n LHS: " + lhs + ".\t RHS: " + rhs);
      System.out.println("X: " + x + "\tY: " + y);
      System.out.println("gcd(" + bigN +", " + (x+y) +": " + factor1 +
                         "\ngcd(" + bigN + ", " + (x-y) + ": " + factor2);
      
      factors.add(Math.abs(factor1));
      factors.add(Math.abs(factor2));
      //check if the factors make other factors!
      if(bigN % factor1 == 0)
        factors.add(Math.abs(bigN/factor1));
      if(bigN % factor2 == 0)
        factors.add(Math.abs(bigN/factor2));
    }
      
    return factors;
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
