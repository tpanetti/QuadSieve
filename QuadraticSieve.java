import java.util.Scanner;
import java.util.ArrayList;

public class QuadraticSieve
{
	private final int SMOOTHNESSBOUND = 1000;
	
	public static void main(String[] args)
	{
	//Grab number N
	int n;
	//take from command line
	if(args.length == 1)
		n = args[0];
	//else take from scanner
	else
	{
		Scanner input = new Scanner(System.in);
		System.out.println("Enter an N: ");
		n = input.nextInt();
	}
	//Find R?
	int R = Math.sqrt(n);	//as long as the sqrt is an integer
	//generate a list of "Big Roots"
	ArrayList<Integer> primes = erathosthenessSieve(R);
	
	
		
	
	
	}
	
	
	public static ArrayList<Integer> eratosthenesSieve(int size) 
	{
		ArrayList<Integer> primes = new ArrayList<Integer>();
		//populate a list with i = 2 to n
		for(int i = 2; i < size; i++) 
			primes.append(i);
		//set p equal to 2?
		int p = 2;
		//loop
		while(Math.Pow(p, 2) < size)
		{
			for(int i = 0; i < primes.size(); i++)
			{
				if(primes.get(i) % p == 0)
					primes.remove(i);
			}
		p++
		}
		return primes;
	}
	
	

}