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
public class QuadraticSieve
{

  private final int SMOOTHNESSBOUND = 1000;
  //DELETE offset (No public vars) 
  private static int offset;
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
    int n;
    int factorbase;
    //take from command line
    if(args.length == 2)
    {
      n = Integer.parseInt(args[0]);
      factorbase = Integer.parseInt(args[1]);
    } 
    //else take from scanner
    else
    {
      Scanner input = new Scanner(System.in);
      System.out.println("Enter an N: ");
      n = input.nextInt();
      factorbase = input.nextInt();
    }
    
    //Find R?
    int R = (int)Math.sqrt(n);  //as long as the sqrt is an integer
    										offset = R; //DELETE PROBABLY

    //generate a list of "Big Roots"
    ArrayList<Integer> primes = eratosthenesSieve(factorbase);
          /* for(int i = 0; i < primes.size(); i++)
            System.out.println(primes.get(i) + " "); */
    //generate a list that satisfies the eqn 
    ArrayList<Integer> bSmooth = findSmoothness(20, n); //MAGIC NUMBER 20
    ArrayList<Integer> residues = calcResiduals(primes, bSmooth);
  
    //Refactor and Gauss
    ArrayList<Pair<ArrayList<Integer>, Integer>> refactoredResidual = refactor(residues, bSmooth, primes);
    //reduce the matrix mod 2
    ArrayList<Pair<ArrayList<Integer>, Integer>> reducedResidual = reduceModTwo(refactoredResidual);
    
    //print reduced Residual matrix
    System.out.println("Reduced matrix is: ");
    int i = 0;
    for(Pair<ArrayList<Integer>,Integer> pair : reducedResidual)
    {
      System.out.print(i + " ");
      System.out.print("[");   
      for(int j : pair.getX())
        System.out.print(j + ", ");
      System.out.println("]");
    i++;
    }
    //perform Gauss-Jordan elimination on the matrix
    ArrayList<ArrayList<Integer>> gaussHits = gauss(reducedResidual);
    //rebuild equation
    HashSet<Integer> factors = rebuildEquation(gaussHits, n);
    System.out.print("The factors are: " );   
    for(int factor : factors)
      System.out.print(factor + ", ");
    
  }
  
  /**
  * This method will perform a Sieve to calculate all primes
  * under the input 
  *
  * @param size This is the size to create the array of primes
  * @return 	An ArrayList of primes under size 
  */        
  public static ArrayList<Integer> eratosthenesSieve(int size) 
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
  
  public static ArrayList<Integer> findSmoothness(int range, int numToFactor)
  {
    //Offset = R. Change later
    ArrayList<Integer> Qs = new ArrayList<Integer>();
    System.out.println("Root is " +offset + "\tindex\t" + "(r+n)^2\t" + "(r+n)^2 - N)");
    for(int i = -range; i <= range; i++) 
    {
      //store log of Q possibly
      int Q = ((int)(Math.pow((offset + i), 2)));
       
      int save = Q - numToFactor;
      Qs.add(save);
      System.out.println(i +"\t" + (Qs.size() - 1) + "\t" + Q + "\t" + save);
    }
  return Qs;
  }
  
  public static ArrayList<Integer> calcResiduals(ArrayList<Integer> primes, ArrayList<Integer> smoothlist)
  {
    //ArrayList<Pair<Integer, Integer>> residuals = new ArrayList<Pair<Integer,Integer>>();

    ArrayList<Integer> copy = new ArrayList<Integer>();
      System.out.println("Size of smoothlist should be 41, size is: " + smoothlist.size());
    for(int i : smoothlist)
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
            int temp = copy.get(j);
            while(temp % i == 0)
            {
              System.out.println("Divide\t" + i + "\t " + (j-20) + "\t" + temp + "\t" + (temp / i));//MAGIC NUMBER
              temp = temp / i;
            }
            copy.set(j, temp);
            j = (j + i) % copy.size();
          } while(j != start);
	}  
      }
      System.out.print("Array [");
      for(int num : copy)
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


  public static ArrayList<Pair<ArrayList<Integer>,Integer>> refactor(ArrayList<Integer> residues, ArrayList<Integer> original,
							ArrayList<Integer> primes)
  {
    //initialize arraylist and fill with zero arryalists
    ArrayList<Pair<ArrayList<Integer>, Integer>> exponents = new ArrayList<Pair<ArrayList<Integer>, Integer>>();

    for(int i = 0; i < residues.size(); i++)
    {
      if(residues.get(i) == 1)
      {
        ArrayList<Integer> exponent = new ArrayList<Integer>();
        for(int index = 0; index < primes.size(); index++)
          exponent.add(0);
	System.out.println("Smooth\t" + (i - 20) + "\t" + i + "\t" + original.get(i)); //MAGIC NUMBER
        int temp = original.get(i);
	if(temp < 0)
	{
	  temp = temp * -1;
	  exponent.set(0,1);
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
        exponents.add(new Pair<ArrayList<Integer>, Integer>(exponent, i));
      }
       
      System.out.print("[");
      for(int whatever : exponents.get(exponents.size() -1).getX())
      {
        System.out.print(whatever +", ");
      }
      System.out.println("]");
    }


    return exponents;
  }
  public static ArrayList<Pair<ArrayList<Integer>, Integer>> reduceModTwo(ArrayList<Pair<ArrayList<Integer>, Integer>> matrix)
  {
   //TODO Print out before and after reduction
   //smooth -1 19 is wrong double check  
    for(int i = 0; i < matrix.size(); i++)
    {
      System.out.print("Row: " + i + "\t[");  
      
      for(int j = 0; j < matrix.get(i).getX().size(); j++)
      {
        System.out.print(matrix.get(i).getX().get(j) + ", ");
        matrix.get(i).getX().set(j, matrix.get(i).getX().get(j) % 2);
  
      }
      System.out.print("]\nReduced\t[");
      for(int temp : matrix.get(i).getX())
        System.out.print(temp + ", ");
      System.out.println("]");
     }
     return matrix;
  }

  
  
  public static ArrayList<ArrayList<Integer>> gauss(ArrayList<Pair<ArrayList<Integer>,Integer>> array)
  {
    System.out.println("Enter Gauss"); 
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
  }
  
  /**
  * This method will rebuild the equation (x
  *
  *
  *
  */
  public static HashSet<Integer> rebuildEquation(ArrayList<ArrayList<Integer>> gaussHits, int bigN)
  {
    int R = (int)(Math.sqrt(bigN));
    HashSet<Integer> factors = new HashSet<Integer>();
    for(ArrayList<Integer> hit : gaussHits)
    {
      int lhs = 1;
      int rhs = 1;
      for(int num : hit)
      {
        lhs *= Math.pow(R + (num - 20), 2);
        rhs *= Math.pow(R + (num - 20), 2) - bigN;
      }
      int x = (int)(Math.sqrt(Math.abs(lhs)));
      int y = (int)(Math.sqrt(Math.abs(rhs)));
      if(x == y)
        continue;
      int factor1 = (int)(gcd(bigN, (x+y)));
      int factor2 = (int)(gcd(bigN, (x-y)));
      System.out.print("Hit\t[");
      for(int number : hit)
      	System.out.print(number + ", "); 
      System.out.println("]\n LHS: " + lhs + ".\t RHS: " + rhs);
      System.out.println("X: " + x + "\tY: " + y);
      System.out.println("gcd(" + bigN +", " + (x+y) +": " + factor1 +
                         "\ngcd(" + bigN + ", " + (x-y) + ": " + factor2);
      factors.add(factor1);
      factors.add(factor2);
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
  public static int gcd(int a, int b)
  {
    while(b != 0)
    {
      int t = b;
      b = a % b;
      a = t;
    }
    return a;
  }



}
