#! /usr/bin/env python
hexcolor = ['#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00','#a65628','#f781bf','#999999']
hexcolor = hexcolor+hexcolor

def qualitative10():
    hexstr = '''#a6cee3
		#1f78b4
		#b2df8a
		#33a02c
		#fb9a99
		#e31a1c
		#fdbf6f
		#ff7f00
		#cab2d6
		#6a3d9a'''
    hexcolor = [s.strip() for s in hexstr.split('\n')]
    return hexcolor

def qualitative8y():
    hexstr = '#e41a1c\n#377eb8\n#4daf4a\n#984ea3\n#ff7f00\n#ffff33\n#a65628\n#f781bf'
    
    hexcolor = [s.strip() for s in hexstr.split('\n')]
    return hexcolor

def qualitative8():
    hexstr = '#e41a1c\n#377eb8\n#4daf4a\n#984ea3\n#ff7f00\n#a65628\n#f781bf\n#999999'
    
    hexcolor = [s.strip() for s in hexstr.split('\n')]
    return hexcolor

def singlehue(color):   
    if color == 'blue':
        hexstr = '''#deebf7
                 #c6dbef
                 #9ecae1
                 #6baed6
                 #4292c6
                 #2171b5
                 #08519c
                 #08306b'''
    elif color == 'green':
        hexstr = '''#e5f5e0
                 #c7e9c0
                 #a1d99b
                 #74c476
                 #41ab5d
                 #238b45
                 #006d2c
                 #00441b'''
    elif color == 'black':
        hexstr = '''#f0f0f0
                 #d9d9d9
                 #bdbdbd
                 #969696
                 #737373
                 #525252
                 #252525
                 #000000'''
    elif color == 'orange':
        hexstr = '''#fee6ce
                 #fdd0a2
                 #fdae6b
                 #fd8d3c
                 #f16913
                 #d94801
                 #a63603
                 #7f2704'''
    elif color == 'purple':
        hexstr = '''#efedf5
                 #dadaeb
                 #bcbddc
                 #9e9ac8
                 #807dba
                 #6a51a3
                 #54278f
                 #3f007d'''
    elif color == 'red':
        hexstr = '''#fee0d2
                 #fcbba1
                 #fc9272
                 #fb6a4a
                 #ef3b2c
                 #cb181d
                 #a50f15
                 #67000d'''
    else:
        print('color "%s" not known.'%color)
        exit(1)
        
    hexcolor = [s.strip() for s in hexstr.split('\n')]
    return hexcolor[::-1]

if __name__ == '__main__':
    print singlehue('red')
    print qualitative8y()
