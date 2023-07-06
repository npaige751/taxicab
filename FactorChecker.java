import java.math.BigInteger;
import java.util.*;
public class FactorChecker {

	static class Factor {
		BigInteger p;
		int c;

		public Factor(BigInteger p) {
			this.p = p;
			this.c = 1;
		}

		public void increment() {
			c ++;
		}

		public BigInteger maxPow() {
			return p.pow(c);
		}

		public String toString() {
			return c > 1 ? p + "^" + c : p.toString();
		}
	}

	public static BigInteger product(ArrayList<Factor> factors, int[] counts) {
		BigInteger res = BigInteger.ONE;
		for (int i=0; i < factors.size(); i++) {
			res = res.multiply(factors.get(i).p.pow(counts[i]));
		}
		return res;
	}

	public static long product(long[] factors, int[] counts) {
		long res = 1;
		for (int i=0; i < factors.length; i++) {
			for (int j=0; j < counts[i]; j++) {
				res *= factors[i];
			}
		}
		return res;
	}

	public static int[] getResidues(int n) {
		boolean[] marks = new boolean[n];
		for (int i=0; i < n; i++) {
			marks[(i*i) % n] = true;
		}
		int c = 0;
		for (int i=0; i < n; i++) {
			if (marks[i]) c++;
		}
		int[] res = new int[c];
		int w = 0;
		for (int i=0; i < n; i++) {
			if (marks[i]) res[w++] = i;
		}
		return res;
	}

	public static final int[] residuesMod256   = getResidues(256);
	public static final int[] residuesMod15015 = getResidues(15015);
	public static final int[] residuesMod12673 = getResidues(12673);
	public static final int[] residuesMod44733 = getResidues(44733);

	public static long check_calls = 0;
	// assumes Y has been premultiplied by 12.
	public static BigInteger[] check(BigInteger X, BigInteger Y) {
		check_calls++;
		BigInteger t = X.multiply(X);
		//t = t.add(t.shiftLeft(1));
		//t = t.add(t).add(t);
		//BigInteger disc = Y.multiply(BigInteger.valueOf(12)).subtract(X.pow(2).multiply(BigInteger.valueOf(3)));
		//BigInteger disc = Y.multiply(BigInteger.valueOf(12)).subtract(t);
		BigInteger disc = Y.subtract(t).subtract(t).subtract(t);
		// disc must be a perfect square
		if (disc.compareTo(BigInteger.ZERO) < 0) return null;
		// check if disc is a quadratic residue mod 256. If evenly distributed, this culls 83% of numbers
		if (Arrays.binarySearch(residuesMod256, disc.intValue() & 0xff) < 0) return null;
		if (Arrays.binarySearch(residuesMod15015, disc.mod(BigInteger.valueOf(15015)).intValue()) < 0) return null;
		if (Arrays.binarySearch(residuesMod12673, disc.mod(BigInteger.valueOf(12673)).intValue()) < 0) return null;
		if (Arrays.binarySearch(residuesMod44733, disc.mod(BigInteger.valueOf(44733)).intValue()) < 0) return null; // this check is about a wash

		BigInteger[] sr = disc.sqrtAndRemainder();
		if (sr[1].equals(BigInteger.ZERO)) {
			BigInteger[] res = new BigInteger[2];
			res[0] = sr[0].divide(BigInteger.valueOf(3)).add(X).divide(BigInteger.valueOf(2));
			res[1] = X.subtract(res[0]);
			if (res[0].compareTo(res[1]) > 0) {
				BigInteger temp = res[0];
				res[0] = res[1];
				res[1] = temp;
			}
			return res;
		}
		return null;
	}

	public static long[] extractSmallprimes(ArrayList<Factor> allFactors) {
		ArrayList<Factor> factors = new ArrayList<>();
		BigInteger spprod = BigInteger.ONE;
		int subsetCount = 1;
		int i = 0;
		while (!allFactors.isEmpty()) {
			BigInteger next = spprod.multiply(allFactors.get(0).maxPow());
			if (next.bitLength() >= 63) break;
			subsetCount *= allFactors.get(0).c + 1;
			spprod = next;
			factors.add(allFactors.remove(0));
		}
		if (spprod.equals(BigInteger.ONE)) return new long[0];
		int maxIndex = factors.size();
		long[] smallprimeChoices = new long[subsetCount];
		long[] smallprimes = new long[maxIndex];
		for (i=0; i < maxIndex; i++) {
			smallprimes[i] = factors.get(i).p.longValue();
		}
		int[] counts = new int[maxIndex];
		int w = 0;
		while (true) {
			smallprimeChoices[w++] = product(smallprimes, counts);
			i = 0;
			counts[0]++;
			while (counts[i] > factors.get(i).c) {
				counts[i] = 0;
				i++;
				if (i == maxIndex) {
					if (w != subsetCount) throw new IllegalStateException("only filled " + w + " of " + subsetCount +  "slots!");
					return smallprimeChoices;
				}
				counts[i]++;
			}
		}
	}

	public static ArrayList<BigInteger[]> check(ArrayList<Factor> factors, BigInteger N) {
		int[] counts = new int[factors.size()];
		ArrayList<BigInteger[]> results = new ArrayList<>();
		while (true) {
			BigInteger X = product(factors, counts);
			BigInteger Y = N.divide(X);
			if (X.compareTo(Y) < 0) {
				BigInteger[] res = check(X, Y);
				if (res != null) results.add(res);
			}
			int i = 0;
			counts[0]++;
			while (counts[i] > factors.get(i).c) {
				counts[i] = 0;
				i++;
				if (i == factors.size()) return results;
				counts[i]++;
			}
		}
	}

	public static ArrayList<BigInteger[]> check2(ArrayList<Factor> factors, BigInteger N) {
		long[] smallprimeChoices = extractSmallprimes(factors);
		long SPP = smallprimeChoices[smallprimeChoices.length-1];	// should be the product of all of the smallprimes
		int[] counts = new int[factors.size()];
		ArrayList<BigInteger[]> results = new ArrayList<>();
		BigInteger ystart = N.divide(BigInteger.valueOf(SPP)).multiply(BigInteger.valueOf(12));
		int sppBits = BigInteger.valueOf(SPP).bitLength();
		int cbrt2n_bits = (N.bitLength() / 3) + 1;
		BigInteger[] spc = new BigInteger[smallprimeChoices.length];
		BigInteger[] spci = new BigInteger[smallprimeChoices.length];
		for (int i=0; i < spc.length; i++) {
			spc[i] = BigInteger.valueOf(smallprimeChoices[i]);
			spci[i] = BigInteger.valueOf(SPP / smallprimeChoices[i]);
		}
		while (true) {
			BigInteger X = product(factors, counts);
			if (X.bitLength() <= cbrt2n_bits) { // X must be less than cbrt(2N)
				BigInteger Y = ystart.divide(X);
				for (int j = 0; j < smallprimeChoices.length; j++) {
					BigInteger x = X.multiply(spc[j]);
					if (x.bitLength() <= cbrt2n_bits) {
						BigInteger y = Y.multiply(spci[j]);
						BigInteger[] res = check(x, y);
						if (res != null) results.add(res);
					}
				}
			}
			if (counts.length == 0) return results;
			int i = 0;
			counts[0]++;
			while (counts[i] > factors.get(i).c) {
				counts[i] = 0;
				i++;
				if (i == factors.size()) return results;
				counts[i]++;
			}
		}
	}

	public static long countSubsets(ArrayList<Factor> factors) {
		long count = 1;
		for (Factor f : factors) {
			count *= f.c+1;
		}
		return count;
	}

	public static void main(String[] args) {
		Scanner sc = new Scanner(System.in);
		int i = 0;
		while (sc.hasNextLine()) {
			String line = sc.nextLine();
			Scanner lsc = new Scanner(line);
			ArrayList<BigInteger> primes = new ArrayList<>();
			BigInteger N = BigInteger.ONE;
			while (lsc.hasNext()) {
				BigInteger p = lsc.nextBigInteger();
				primes.add(p);
				N = N.multiply(p);
			}
			Collections.sort(primes);
			ArrayList<Factor> factors = new ArrayList<>();
			BigInteger prev = null;
			Factor prevFac = null;
			for (BigInteger p : primes) {
				if (p.equals(prev)) {
					prevFac.increment();
				} else {
					Factor f = new Factor(p);
					factors.add(f);
					prev = p;
					prevFac = f;
				}
			}
			System.out.print("\rChecking i = " + i + "; N = " + N + " (" + countSubsets(factors) + " subsets)");
			ArrayList<BigInteger[]> results = check2(factors, N);
			if (results.size() != 2) {
				System.out.print("\rindex " + i + "; " + N + " is " + results.size() + "-way: N");
				for (BigInteger[] res : results) {
					System.out.print(" = " + res[0] + "^3 + " + res[1] + "^3");
				}
				System.out.println();
			}
			i++;
		}
		System.out.println("\nCheck calls: " + check_calls);
	}
}
