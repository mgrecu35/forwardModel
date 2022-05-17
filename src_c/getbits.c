/* getbits:  get n bits from position p */
unsigned long getbits(unsigned long x, int p, int n)
{
     return (x >> (p+1-n)) & ~(~0 << n);
}
