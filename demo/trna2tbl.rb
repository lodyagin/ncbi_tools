#! /usr/bin/ruby

# subroutines

def firstnum str
  a = str.index('(')
  b = str.index('-')
  loc = str.slice(a + 1 .. b - 1)
  loc.to_i
end

def lastnum str
  c = str.index('-')
  d = str.index(')')
  loc = str.slice(c + 1 .. d - 1)
  loc.to_i
end

def getaa str
  x = str.index(' ')
  y = str.index("\t")
  loc = str.slice(x + 1 .. y - 1)
  loc.to_s
end

def getid str
  z = str.index('.trna')
  loc = str.slice(0 .. z - 1)
  loc.to_s
end

# main program

idnotsent = true

st = 0;
sp = 0;
inst = 0
insp = 0
acst = 0
acsp = 0
id = ''
aa = ''
pseudo = false

name = gets

while name
  name = name.chomp

  if name.index('Seq:') != nil

    if idnotsent
      puts '>Features ' + id + ' tRNAscan-SE'
      idnotsent = false
    end

    if inst == 0 and insp == 0
      puts st.to_s + "\t" + sp.to_s + "\t" + "tRNA"
    else
      inst = inst - 1
      insp = insp + 1
      puts st.to_s + "\t" + inst.to_s + "\t" + "tRNA"
      puts insp.to_s + "\t" + sp.to_s
    end

    puts "\t\t\tproduct\t" + aa

    if pseudo
      puts "\t\t\tpseudo"
    else
      if acst != 0 and acsp != 0
        if acst < acsp
          puts "\t\t\tanticodon\t(pos:" + acst.to_s + ".." + acsp.to_s + ",aa:" + aa + ")"
        else
          puts "\t\t\tanticodon\t(pos:complement(" + acsp.to_s + ".." + acst.to_s + "),aa:" + aa + ")"
        end
      end
    end

    puts "\t\t\tgene\t-"

    st = 0;
    sp = 0;
    inst = 0
    insp = 0
    acst = 0
    acsp = 0
    id = ''
    aa = ''
    pseudo = false

  elsif name.index('Type:') != nil

    aa = getaa name
    if aa == 'Undet' or aa == 'Sup'
      aa = 'OTHER'
    end

    loc = name.slice(/[(].*[)]/)
    acst = firstnum loc
    acsp = lastnum loc

  elsif name.index('Length') != nil

    id = getid name

    loc = name.slice(/[(].*[)]/)
    st = firstnum loc
    sp = lastnum loc

  elsif name.index('Possible intron:') != nil

    loc = name.slice(/[(].*[)]/)
    inst = firstnum loc
    insp = lastnum loc

  elsif name.index('Possible pseudogene:') != nil

    pseudo = true

  end
  name = gets
end

