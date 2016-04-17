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
		//take from command line
		if(args.length == 1)
			n = Integer.parseInt(args[0]);
		//else take from scanner
		else
		{
			Scanner input = new Scanner(System.in);
			System.out.println("Enter an N: ");
			n = input.nextInt();
		}
		//Find R?
		int R = (int)Math.sqrt(n);	//as long as the sqrt is an integer
		//generate a list of "Big Roots"
		ArrayList<Integer> primes = eratosthenesSieve(n);
					/* for(int i = 0; i < primes.size(); i++)
						System.out.println(primes.get(i) + " "); */
		//generate a list that satisfies 
		ArrayList<Integer> bSmooth = findSmoothness(R, n, 100);
		//ArrayList<Tuple<Integer, Integer>> residues = calcResidue(primes, bSmooth);
	
	
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
	
	public static ArrayList<Integer> findSmoothness(int floorN, int factor, int n)
	{
		ArrayList<Integer> Qs = new ArrayList<Integer>();
		for(int i = -n; i <= n; i++) 
		{
			//store log of Q possibly
			int Q = ((int)(Math.pow((floorN + n), 2))) - factor;
			Qs.add(Q);
		}
	return Qs;
	}
	
	//public static ArrayList<Tuple<Integer, Integer>> calcResidue(

}
