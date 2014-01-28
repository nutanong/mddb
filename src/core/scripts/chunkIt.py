

def fixed_num(seq, num):
  avg = len(seq) / float(num)
  out = []
  last = 0.0
  while last < len(seq):
    out.append(seq[int(last):int(last + avg)])
    last += avg
  return out

def fixed_size(seq, size):
    for i in xrange(0, len(seq), size):
        yield seq[i:i+size]

