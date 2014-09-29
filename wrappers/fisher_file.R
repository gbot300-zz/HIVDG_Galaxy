library('stats')
my.fisher = function(x,y){
	fisher.test(array(c(x,y),dim=c(length(x),2)),workspace = 2e8)
}
