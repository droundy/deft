set -ev

#python figs/yaml-comparison.py /home/jordan/sad-monte-carlo/ data/samc-1e5-256-reference-lndos.dat s000/periodic-ww1.50-ff0.17-N256 256 915 509 true samc-1e4-256 wl-256 samc-1e7-256 samc-1e8-256
#python figs/yaml-comparison.py /home/jordan/sad-monte-carlo/ data/samc-1e5-256-reference-lndos.dat s000/periodic-ww1.50-ff0.17-N256 256 737 509 true samc-1e3-256

#python figs/yaml-comparison.py /home/droundy/src/sad-monte-carlo/ data/samc-1e5-256-reference-lndos.dat s000/periodic-ww1.50-ff0.17-N256 256 915 509 true samc-1e5-256 samc-1e6-256 sad-256
#python figs/yaml-comparison.py /home/jordan/sad-monte-carlo/ data/samc-1e5-256-reference-lndos.dat s000/periodic-ww1.50-ff0.17-N256 256 915 509 true wl-inv-t-256

#python figs/yaml-comparison.py /home/droundy/src/sad-monte-carlo/ data/samc-1e3-T13-reference-lndos.dat s000/periodic-ww1.30-ff0.30-N50 50 248 120 true sad-T13 sad-slow-T13

#python figs/yaml-comparison.py /home/jordan/sad-monte-carlo/ data/samc-1e3-T13-reference-lndos.dat s000/periodic-ww1.30-ff0.30-N50 50 242 120 true samc-1e3-50
#python figs/yaml-comparison.py /home/jordan/sad-monte-carlo/ data/samc-1e3-T13-reference-lndos.dat s000/periodic-ww1.30-ff0.30-N50 50 248 120 true samc-1e4-50 samc-1e5-50 samc-1e6-50 samc-1e7-50 wl-50

#python figs/yaml-comparison.py /home/jordan/sad-monte-carlo/ data/samc-1e3-T13-reference-lndos.dat s000/periodic-ww1.30-ff0.30-N50 50 248 120 true samc-1e3-50-slow
#python figs/yaml-comparison.py /home/jordan/sad-monte-carlo/ data/samc-1e3-T13-reference-lndos.dat s000/periodic-ww1.30-ff0.30-N50 50 248 146 true samc-1e4-50-slow
#python figs/yaml-comparison.py /home/jordan/sad-monte-carlo/ data/samc-1e3-T13-reference-lndos.dat s000/periodic-ww1.30-ff0.30-N50 50 248 120 true samc-1e5-50-slow samc-1e6-50-slow samc-1e7-50-slow wl-50-slow

#python figs/yaml-comparison.py /home/jordan/sad-monte-carlo/ data/samc-1e3-T13-reference-lndos.dat s000/periodic-ww1.30-ff0.30-N50 50 248 120 true wl-inv-t-50 wl-inv-t-50-slow

#python figs/yaml-comparison.py /home/jordan/sad-monte-carlo/ data/samc-1e3-T13-reference-lndos.dat s000/periodic-ww1.30-ff0.30-N50 50 248 120 true sa-50 sa-50-slow
#python figs/yaml-comparison.py /home/jordan/sad-monte-carlo/ data/samc-1e5-256-reference-lndos.dat s000/periodic-ww1.50-ff0.17-N256 256 915 509 true sa-256


python figs/parse-yaml-out.py /home/jordan/weniger-cp-data/sad-50-slow-s1 50 sad-50-slow
python figs/parse-yaml-out.py /home/jordan/weniger-cp-data/sad-50-s1 50 sad-50
python figs/parse-yaml-out.py /home/jordan/weniger-cp-data/wl-50-s1 50 wl-50
python figs/parse-yaml-out.py /home/jordan/weniger-cp-data/wl-inv-t-50-s1 50 wl-inv-t-50

python figs/parse-yaml-out.py /home/jordan/weniger-cp-data/sad-256-s1 256 sad-256
python figs/parse-yaml-out.py /home/jordan/weniger-cp-data/wl-256-s1 256 wl-50
python figs/parse-yaml-out.py /home/jordan/weniger-cp-data/wl-inv-t-256-s1 256 wl-inv-t-256

python figs/yaml-multi-comparison.py /home/jordan/weniger-cp-data/ data/samc-1e3-T13-reference-lndos.dat s000/periodic-ww1.30-ff0.30-N50 50 248 120 true 4 sad-50 wl-50 wl-inv-t-50 samc-1e3-50 samc-1e4-50 samc-1e5-50 samc-1e6-50 samc-1e7-50
python figs/yaml-multi-comparison.py /home/jordan/weniger-cp-data/ data/samc-1e3-T13-reference-lndos.dat s000/periodic-ww1.30-ff0.30-N50 50 248 120 true 4 sad-50-slow wl-50-slow wl-inv-t-50-slow samc-1e3-50-slow samc-1e4-50-slow samc-1e5-50-slow samc-1e6-50-slow samc-1e7-50-slow
python figs/yaml-multi-comparison.py /home/jordan/weniger-cp-data/ data/samc-1e5-256-reference-lndos.dat s000/periodic-ww1.50-ff0.17-N256 256 915 509 true 4 sad-256 wl-256 wl-inv-t-256 samc-1e3-256 samc-1e4-256 samc-1e5-256 samc-1e6-256 samc-1e7-256 samc-1e8-256
