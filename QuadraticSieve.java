import java.util.Scanner;
import java.util.ArrayList;
/*
**This program will run a quadratic sieve to find
**some numbers?
**@author Thomas Panetti & Noemi
**
**
*/
public class QuadraticSieve
{
  //TODO: import proper tuple

  private final int SMOOTHNESSBOUND = 1000;
  private static int offset;
  private class Pair<X, Y> 
  {
    public final X x;
    public final Y y;
    public Pair(X x, Y y) 
    {
      this.x = x;
      this.y = y;
    }
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
    ArrayList<Integer> bSmooth = findSmoothness(R, n);
    ArrayList<Integer> residues = calcResiduals(primes, bSmooth);
  
    //Refactor and Gauss
    ArrayList<Integer> refactoredResidual = refactor(residues, bSmooth); 
  }
  
          
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
    return primes;
  }
  
  public static ArrayList<Integer> findSmoothness(int range, int numToFactor)
  {
    ArrayList<Integer> Qs = new ArrayList<Integer>();
    System.out.println("Root is " +range + "\tindex\t" + "(r+n)^2\t" + "(r+n)^2 - N)");
    for(int i = -range; i <= range; i++) 
    {
      //store log of Q possibly
      int Q = ((int)(Math.pow((range + i), 2)));
       
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
    for(int i : smoothlist)
      copy.add(i);
    int start;
    for(int i : primes)
    {
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
              System.out.println("Divide\t" + i + "\t " + (j-offset) + "\t" + temp + "\t" + (temp / i));
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


  public static ArrayList<ArrayList<<Integer>> refactor(ArrayList<Integer> residues, ArrayList<Integer> original)
  {
    return null;
  }

  public static ArrayList<ArrayList<Integer>> Gauss(ArrayList<Intger
  {
    return null;
  }


}
