package io.jr576.utils;

import java.util.List;
import java.util.stream.IntStream;
import java.util.stream.Stream;

public class FFT {
    final long prime;
    final long primRoot;
    final int order;
    final long[] rootPowers;
    private final int recursiveParallelBreakpoint = 1 << 6;
    private final int useRecursiveAlgorithm = 1 << 10;
    public static final List<Integer> goodPrimes = List.of(469762049, 167772161, 754974721, 377487361, 595591169, 645922817, 880803841, 897581057, 998244353);

    private long[] initializeRootPowers() {
        long[] rootPowers = new long[order];
        rootPowers[0] = 1;
        for (int i = 1; i < order; i++)
            rootPowers[i] = rootPowers[i - 1] * primRoot % prime;
        return rootPowers;
    }

    public FFT(int prime, int primRoot, int order) {
        if (Integer.highestOneBit(order) != order)
            throw new IllegalArgumentException("Order must be power of two.");
        this.prime = prime;
        this.primRoot = primRoot;
        this.order = largerTwoPow(order);
        this.rootPowers = initializeRootPowers();
    }

    public FFT(int prime) {
        this.prime = prime;
        int maxOrder = prime - 1;
        int maxEvenOrder = Integer.lowestOneBit(maxOrder);
        int oddOrder = maxOrder / maxEvenOrder;
        long primRoot = 2;
        while (pow(primRoot, maxOrder / 2) == 1)
            primRoot++;
        primRoot = pow(primRoot, oddOrder);
        this.order = maxEvenOrder;
        this.primRoot = primRoot;
        this.rootPowers = initializeRootPowers();
    }

    public FFT(int prime, int desiredOrder) {
        this.prime = prime;
        this.order = desiredOrder;
        desiredOrder = largerTwoPow(desiredOrder);
        int maxOrder = prime - 1;
        if (maxOrder < desiredOrder)
            throw new IllegalArgumentException("The prime doesn't support order this large.");
        long primRoot = 2;
        while (pow(primRoot, maxOrder / 2) == 1)
            primRoot++;
        primRoot = pow(primRoot, maxOrder / desiredOrder);
        this.primRoot = primRoot;
        this.rootPowers = initializeRootPowers();
    }

    private long pow(long base, long exponent) {
        long result = 1;
        while (exponent > 0) {
            if ((exponent & 1) == 1)
                result = result * base % prime;
            base = base * base % prime;
            exponent >>= 1;
        }
        return result;
    }

    private long[] extendAndCopy(long[] arr, int size) {
        long[] newArr = new long[size];
        System.arraycopy(arr, 0, newArr, 0, arr.length);
        return newArr;
    }

    //smallest power of two >= n
    private int largerTwoPow(int n) {
        if (Integer.bitCount(n) == 1)
            return n;
        else {
            int highestOneBit = Integer.highestOneBit(n);
            return highestOneBit << 1;
        }
    }

    private long[] fftRecursive(long[] pol, int size, int rootIndex) {
        if (size > order)
            throw new IllegalArgumentException("Size too large.");
        if (size == 1)
            return pol;
        long[] even = new long[size / 2];
        long[] odd = new long[size / 2];
        for (int i = 0; i < size / 2; i++) {
            even[i] = pol[2 * i];
            odd[i] = pol[2 * i + 1];
        }
        if (size > recursiveParallelBreakpoint)
            Stream.of(even, odd)
                  .parallel()
                  .forEach(array -> fftRecursive(array, size / 2, 2 * rootIndex));
        else {
            fftRecursive(even, size / 2, 2 * rootIndex);
            fftRecursive(odd, size / 2, 2 * rootIndex);
        }
        for (int k = 0; k < size / 2; k++) {
            long evenPart = even[k];
            long oddPart = rootPowers[k * rootIndex] * odd[k] % prime;
            pol[k] = (evenPart + oddPart) % prime;
            pol[k + size / 2] = (evenPart + prime - oddPart) % prime;
        }
        return pol;
    }

    private long[] fftIterative(long[] pol, int size) {
        if (size > order)
            throw new IllegalArgumentException("Size too large.");
        if (size == 1)
            return pol;
        bitReverse(pol);
        for (int twoPow = 2; twoPow <= size; twoPow *= 2) {
            int twoPowFinal = twoPow;
            int rootIndexIncrement = rootPowers.length / twoPow;
            IntStream.range(0, size / twoPowFinal).parallel().forEach(i -> {
                int k = i * twoPowFinal;
                for (int j = 0, rootIndex = 0; j < twoPowFinal / 2; j++, rootIndex += rootIndexIncrement) {
                    long root = rootPowers[rootIndex];
                    long evenPart = pol[k + j];
                    long oddPart = root * pol[k + j + twoPowFinal / 2] % prime;
                    pol[k + j] = (evenPart + oddPart) % prime;
                    pol[k + j + twoPowFinal / 2] = (evenPart + prime - oddPart) % prime;
                }
            });
        }
        return pol;
    }

    private long[] bitReverse(long[] pol) {
        int size = pol.length;
        int shift = 1 + Integer.numberOfLeadingZeros(size);
        IntStream.range(0, size).parallel().forEach(i -> { //paralell??
            int j = Integer.reverse(i) >>> shift;
            if (j > i) {
                long temp = pol[i];
                pol[i] = pol[j];
                pol[j] = temp;
            }
        });
        return pol;
    }

    private long[] fft(long[] pol, int size) {
        long[] polToTransform = extendAndCopy(pol, size);
        if (size < useRecursiveAlgorithm)
            return fftRecursive(polToTransform, size, order / size);
        else
            return fftIterative(polToTransform, size);
    }

    public long[] fft(long[] pol) {
        int size = largerTwoPow(pol.length);
        return fft(pol, size);
    }

    private long[] fftInverse(long[] pol) {
        long[] polTransform = fft(pol);
        long[] fftInverse = new long[polTransform.length];
        long sizeInverse = pow(polTransform.length, prime - 2);
        fftInverse[0] = polTransform[0] * sizeInverse % prime;
        for (int i = 1; i < polTransform.length; i++)
            fftInverse[i] = polTransform[polTransform.length - i] * sizeInverse % prime;
        return fftInverse;
    }

    public long[] multiply(long[] polA, long[] polB) {
        int size = largerTwoPow(polA.length + polB.length - 1);
        long[] polATransformed = fft(extendAndCopy(polA, size), size);
        long[] polBTransformed = fft(extendAndCopy(polB, size), size);
        long[] product = new long[size];
        for (int i = 0; i < size; i++)
            product[i] = polATransformed[i] * polBTransformed[i] % prime;
        return fftInverse(product);
    }

    public long[] multiplyCutoff(long[] polA, long[] polB, int cutoff) {
        long[] product = multiply(polA, polB);
        long[] productCutoff = new long[cutoff + 1];
        System.arraycopy(product, 0, productCutoff, 0, cutoff + 1);
        return productCutoff;
    }

    public long[] square(long[] pol) {
        int size = largerTwoPow(2 * pol.length - 1);
        long[] polTransformed = fft(pol, size);
        long[] square = new long[size];
        for (int i = 0; i < size; i++)
            square[i] = polTransformed[i] * polTransformed[i] % prime;
        return fftInverse(square);
    }

    public long[] squareCutoff(long[] pol, int cutoff) {
        long[] square = square(pol);
        if (square.length == cutoff + 1)
            return square;
        long[] squareCutoff = new long[cutoff + 1];
        System.arraycopy(square, 0, squareCutoff, 0, cutoff + 1);
        return squareCutoff;
    }

    public long[] pow(long[] pol, int exponent) {
        int size = largerTwoPow(exponent * (pol.length - 1) + 1);
        long[] polTransformed = fft(pol, size);
        long[] power = new long[size];
        for (int i = 0; i < size; i++)
            power[i] = pow(polTransformed[i], exponent);
        return fftInverse(power);
    }


    public long[] powCutoff(long[] base, long exponent, int cutoff) {
        long[] repeatedSquare = base;
        while (exponent % 2 == 0) {
            repeatedSquare = squareCutoff(repeatedSquare, cutoff);
            exponent /= 2;
        }
        long[] power = repeatedSquare.clone();
        repeatedSquare = squareCutoff(repeatedSquare, cutoff);
        exponent /= 2;
        while (exponent > 0) {
            if (exponent % 2 == 1)
                power = multiplyCutoff(power, repeatedSquare, cutoff);
            repeatedSquare = squareCutoff(repeatedSquare, cutoff);
            exponent /= 2;
        }
        return power;
    }

    private long[] transformSlow(long[] pol, int size, int rootIndex) {
        long[] eval = new long[size];
        for (int j = 0; j < size; j++) {
            long result = 0;
            long root = rootPowers[j * rootIndex];
            for (int i = pol.length - 1; i >= 0; i--) {
                result = (result * root + pol[i]) % prime;
            }
            eval[j] = result;
        }
        return eval;
    }

}
