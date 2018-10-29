
<?php




 /*
       需求：求解函数 f(x) = x + 10*sin(5*x) + 7*cos(4*x) 在区间[0,9]的最大值。
   
   
       以我们的目标函数 f(x) = x + 10sin(5x) + 7cos(4x), x∈[0,9] 为例。
       假如设定求解的精度为小数点后4位，可以将x的解空间划分为 (9-0)×(1e+4)=90000个等分。
       2^16<90000<2^17，需要17位二进制数来表示这些解。换句话说，一个解的编码就是一个17位的二进制串。
       一开始，这些二进制串是随机生成的。
      一个这样的二进制串代表一条染色体串，这里染色体串的长度为17。
      对于任何一条这样的染色体chromosome，如何将它复原(解码)到[0,9]这个区间中的数值呢？
  
      对于本问题，我们可以采用以下公式来解码：x = 0 + decimal(chromosome)×(9-0)/(2^17-1)
  
  */
  header("Content-Type:text/html;CharSet=UTF-8");
  //初始化参数
  $elitism = true;    //是否精英选择
  $population_size = 100; //种群大小
  $chromosome_size = 17;  //染色体长度
  $generation_size = 200; //最大迭代次数
  $cross_rate = 0.6;      //交叉概率
  $mutate_rate = 0.01;    //变异概率
  
  
  //初始化对象，执行算法
  $ga_object = new GAEngine($population_size, $chromosome_size, $generation_size, $cross_rate, $mutate_rate, $elitism);
  
  $res = $ga_object->run();
  
  
  //打印结果
  echo "迭代{$generation_size}代，最佳个体如下<br>";
  echo "染色体：{$res['m']}<br>";
  echo "对应的X：{$res['q']}<br>";
  echo "结果：{$res['n']}<br>";
  echo "最佳个体出现的代：【{$res['p']}】<br>";
  
  
  
  
  
  //遗传算法主要部分
  class GAEngine
  {
      //初始化参数
      public $elitism;            //是否精英选择
      public $population_size;    //种群大小
      public $chromosome_size;    //染色体长度
      public $generation_size;    //最大迭代次数
      public $cross_rate;         //交叉概率
      public $mutate_rate;        //变异概率
  
      //全局参数
      public $G;  //当前迭代次数
      public $fitness_value;  //当前代适应度矩阵
      public $best_fitness;   //历代最佳适应值
      public $fitness_sum;    //前i个个体的适应度之和
      public $fitness_average;    //历代平均适应值矩阵
      public $best_individual;    //历代最佳个体
      public $best_generation;    //最佳个体出现代
      public $upper_bound = 9;    //自变量的区间上限
      public $lower_bound = 0;    //自变量的区间下限
  
      //对象
      public $population_obj;     //种群对象
  
  
      public $populations;    //种群数组
 
  
      //初始化
      public function __construct($population_size, $chromosome_size, $generation_size, $cross_rate, $mutate_rate, $elitism)
      {
          $this->elitism = $elitism;
          $this->population_size = $population_size;
          $this->chromosome_size = $chromosome_size;
          $this->generation_size = $generation_size;
          $this->cross_rate = $cross_rate;
          $this->mutate_rate = $mutate_rate;
  
          //初始化种群
          $this->population_obj = new Populations($population_size, $chromosome_size);
  		
      }
  
      //【执行,返回结果】
      public function run(){
          //初始化种群
          $this->populations = $this->population_obj->initPop();
  		var_dump( $this->populations);
          //开始迭代
          for($i = 1; $i < $this->generation_size; $i ++){
              $this->G = $i;
              $this->fitness();  //计算适应度
              $this->rank();     //对个体按适应度大小进行排序
              $this->selection();//选择操作
              $this->crossover();    //交叉操作
             $this->mutation();     //变异操作
         }
 
 
         //最后，求出二进制基因对应的实数x，以及求出实数对应的结果
         $x = $this->upper_bound - bindec($this->best_individual)*($this->upper_bound - $this->lower_bound)/(pow(2, $this->chromosome_size) - 1);
 
         return [
             'm' => $this->best_individual,  //最佳个体
             'n' => $this->best_fitness,     //最佳适应度
             'p' => $this->best_generation,  //最佳个体出现的代
             'q' => $x                       //最佳个体基因对应的实数
         ];
 
     }
 
     //【计算适应度】
     public function fitness(){
         //所有个体适应度初始化为0
         //遍历每个个体的基因，求出对应的实数，然后再计算出结果，即适应度
         for($i = 0; $i < $this->population_size; $i ++){
 //            $gens = strrev($this->populations[$i]['gens']);//染色体字符串与实际的自变量x二进制串顺序是相反的
             $x = $this->upper_bound - bindec($this->populations[$i]['gens'])*($this->upper_bound - $this->lower_bound)/(pow(2, $this->chromosome_size) - 1);
 
             $fx = $x + 10*sin(5*$x) + 7*cos(4*$x);//函数值
 
             $this->fitness_value[$i] = $fx;
 
         }
     }
 
     //【排序】
     //对个体按适应度的大小进行排序,并且保存最佳个体
     public function rank(){
         //$this->fitness_value[] 保存适应度的数组，根据此数组进行排序，还要修改对应个体的位置
         //冒泡，从小到大排序
         for($i = 0; $i < $this->population_size -1; $i ++){
             for($j = $i + 1; $j < $this->population_size; $j ++){
                 if($this->fitness_value[$i] > $this->fitness_value[$j]){
                     //交换适应度
                     $tmp = $this->fitness_value[$i];
                     $this->fitness_value[$i] = $this->fitness_value[$j];
                     $this->fitness_value[$j] = $tmp;
 
                     //交换种群个体位置
                     $tmp = $this->populations[$i];
                     $this->populations[$i] = $this->populations[$j];
                     $this->populations[$j] = $tmp;
                 }
             }
         }
 
         //计算前i个给个体的适应度之和
         $fitness_sum[0] = $this->fitness_value[0];
         for($i = 1; $i < $this->population_size; $i ++){
             $fitness_sum[$i] = $fitness_sum[$i - 1] +$this->fitness_value[$i];
         }
         $this->fitness_sum = $fitness_sum;
 
         //第G代迭代， 个体的平均适应度
         $this->fitness_average[$this->G] = ($fitness_sum[$this->population_size - 1] / $this->population_size);
 
         //更新最大适应度和对应的迭代次数，保存最佳个体（最佳个体适应度最大）
         if($this->fitness_value[$this->population_size - 1] > $this->best_fitness){
             $this->best_fitness = $this->fitness_value[$this->population_size - 1];
             $this->best_generation = $this->G;
             $this->best_individual = $this->populations[$this->population_size -1]['gens'];
         }
 
     }
 
     //【选择】
     public function selection(){
         $population_new = [];//保存被选中的个体
 
         //二分查找实现轮盘赌功能
         for($i = 0; $i < $this->population_size; $i ++){
             $r = (rand(0,10)/10 ) * $this->fitness_sum[$this->population_size - 1];//生成一个随机数，在[0,总适应度] 之间
             $first = 1;
             $last = $this->population_size;
             $mid = round( ($last + $first) / 2 );
             $idx = -1;
 
 
             //排中法选择个体
             while(($first <= $last) && ($idx == -1)){
                 if($r > $this->fitness_sum[$mid]){
                     $first = $mid;
                 }elseif ( $r < $this->fitness_sum[$mid]){
                     $last = $mid;
                 }else{
                     $idx = $mid;
                     break;
                 }
                 $mid = round( ($last + $first) / 2 );
                 if(($last - $first) == 1){
                     $idx = $last;
                     break;
                 }
 
             }
 
             //产生新的个体
             //echo $idx.'=';
             $population_new[$i]['gens'] = $this->populations[$idx]['gens'];
         }
 
         //是否精英选择
         if($this->elitism){
             $p = $this->population_size - 1;
         }else{
             $p = $this->population_size;
         }
 
         for($i = 0; $i < $p; $i ++){
             if(isset($population_new[$i]))
                 $this->populations[$i] = $population_new[$i];
         }
 
     }
 
    //【交叉】
     public function crossover(){
         //步长为2, 遍历种群
         for($i = 0; $i < $this->population_size; $i+=2){
             //rand < 交叉概率，对2个个体的染色体进行交叉操作
             $r = $this->randFloat();//产生0~1的随机数
             if($r < $this->cross_rate){
 //                $r =$this->randFloat();
 //                $cross_pos = round($r * $this->chromosome_size);
                 $cross_pos = rand(0, $this->chromosome_size);
                 if($cross_pos ==0 || $cross_pos == $this->chromosome_size)
                     continue;
                 //对 cross_pos及之后的二进制串进行交换
                 $x = $this->populations[$i]['gens'];
                 $y =  $this->populations[$i+1]['gens'];
                 $tmp1 = substr($x,0, $cross_pos).substr($y,$cross_pos);
                $tmp2 = substr($y,0,$cross_pos).substr($x, $cross_pos);
 
                 $this->populations[$i]['gens'] = $tmp1;
                 $this->populations[$i+1]['gens'] = $tmp2;
             }
         }
     }
 
     //【变异】
     public function mutation(){
         for($i =0; $i < $this->population_size; $i ++){
            if($this->randFloat() < $this->mutate_rate){
                $mutate_pos = rand(0,$this->chromosome_size -1);
                $this->populations[$i]['gens'][$mutate_pos] = 1 - $this->populations[$i]['gens'][$mutate_pos];
            }
 
         }
     }
 
  
     public function show($data){
         echo '<pre>';
         var_dump($data);
         echo '<hr>';
     }
 
     //随机产生0-1的小数
     function randFloat($min=0, $max=1){
         return $min + mt_rand()/mt_getrandmax() * ($max-$min);
     }
 
 }
 
 //种群
 class Populations
 {
     public $population_size;//种群大小
     public $chromosome_size;//染色体长度
 
     //初始化参数
     public function __construct($population_size, $chromosome_size)
     {
         $this->population_size = $population_size;
         $this->chromosome_size = $chromosome_size;
     }
 
     //初始化种群
     public function initPop(){
         $pop = [];
         for($i = 0; $i < $this->population_size; $i ++){
         	$pop[$i] = [];
         	$pop[$i]['gens'] ='';
             for($j = 0; $j < $this->chromosome_size; $j ++){
                 $pop[$i]['gens'] .=  rand(0, 10) % 2;//产生1-0随机数
             }
         }
 
         return $pop;
     }
 
 }


// $res = new Populations(100,17);

//  print_r($res->initPop());













// $url = "http://test.yusj.vip/?s=Video.drama.test&token=NhiDYwIQB7K3DXlh";


// $data = json_encode(['a'=>'b','c'=>'d']);
// echo date('Y-m-d',strtotime('-1 day'));
// $id = '4200000146201806145893735731';
// $id=str_replace(",","",number_format($id));
// echo $id;
// $id = substr($id,10,4);
// echo $id;
// $create_time = substr($id,10,4).'-'.substr($id,14,2).'-'.substr($id,16,2);
//  echo $create_time;
// echo round(3.256,2);

// $url = "http://test.yusj.vip/?s=Order.receive.notifycurrencys";

// $data = [
// 'data'=>[
//  'data'=>base64_encode(json_encode([
//      'out_trade_no'=>2018092779049,
//      'trade_no'=>123456789,
//      'type'=>'wxwb',
//      'desc'=>'支付成功'
//      ])),
//  'pay_info'=>base64_encode(json_encode($_POST)),
//  'sign'=>'798bdd1ec31d173a'
// ]
// ];

// $output = curl_request($url,$data);

// var_dump($output);

// //curl调用
// function curl_request($url,$post='',$cookie='', $returnCookie=0){
//           $curl = curl_init();
//           curl_setopt($curl, CURLOPT_URL, $url);
//           curl_setopt($curl, CURLOPT_USERAGENT, 'Mozilla/5.0 (compatible; MSIE 10.0; Windows NT 6.1; Trident/6.0)');
//           curl_setopt($curl, CURLOPT_FOLLOWLOCATION, 1);
//           curl_setopt($curl, CURLOPT_AUTOREFERER, 1);

//           curl_setopt($curl, CURLOPT_REFERER, "http://XXX");
//           if($post) {
//              curl_setopt($curl, CURLOPT_POST, 1);
//              curl_setopt($curl, CURLOPT_POSTFIELDS, http_build_query($post));
//          }
//          if($cookie) {
//              curl_setopt($curl, CURLOPT_COOKIE, $cookie);
//          } 
//          curl_setopt($curl, CURLOPT_HEADER, $returnCookie);
//          curl_setopt($curl, CURLOPT_TIMEOUT, 10);
//         curl_setopt($curl, CURLOPT_RETURNTRANSFER, 1);
//          $data = curl_exec($curl);
//          if (curl_errno($curl)) {
//              return curl_error($curl);
//          }
//          curl_close($curl);
//          if($returnCookie){
//              list($header, $body) = explode("\r\n\r\n", $data, 2);
//              preg_match_all("/Set\-Cookie:([^;]*);/", $header, $matches);
//              $info['cookie']  = substr($matches[1][0], 1);
//              $info['content'] = $body;
//              return $info;
//          }else{
//             return $data;
//          }
//  }




// $a='\<xml\>\<';
// // echo $a;                  
// var_dump($a);

//  $a = http_request($url,$data);
// echo '<pre>';
//  $b = json_decode($a,true);
//  echo $b['ret'];
// echo '</pre>';
//  function http_request( $url, $post = '', $timeout = 5 ){ 
// if( empty( $url ) ){
//  return ;
// }
// $ch = curl_init();
// curl_setopt($ch, CURLOPT_URL, $url);
// //头部文件
// curl_setopt($ch, CURLOPT_HEADER, 0);
// //数据流输出
// curl_setopt($ch, CURLOPT_RETURNTRANSFER, 1);

// //取消ssl验证
// curl_setopt($ch, CURLOPT_SSL_VERIFYPEER, FALSE);
// curl_setopt($ch, CURLOPT_SSL_VERIFYHOST, FALSE);
 
// if( $post != '' && !empty( $post ) ){
// 	//post方式提交
//  curl_setopt($ch, CURLOPT_POST, 1);
//  //提交数据
//  curl_setopt($ch, CURLOPT_POSTFIELDS, $post);
//  //数据格式
//  curl_setopt($ch, CURLOPT_HTTPHEADER, array('Content-Type: application/json', 'Content-Length: ' . strlen($post)));
// }
// //提交数据
// curl_setopt($ch, CURLOPT_TIMEOUT, $timeout);
// //抓取返回
// $result = curl_exec($ch);
// curl_close($ch);
// return $result;
// }

//ip地址对应的指向地址
// function getIPLoc_sina($ip){

//    $url = 'http://int.dpool.sina.com.cn/iplookup/iplookup.php?format=json&ip='.$queryIP;
//
//    $ch = curl_init($url);
//
//    //curl_setopt($ch,CURLOPT_ENCODING ,'utf8');
//
//    curl_setopt($ch, CURLOPT_TIMEOUT, 10);
//
//    curl_setopt($ch, CURLOPT_RETURNTRANSFER, true) ; // 获取数据返回
//
//    $location = curl_exec($ch);
//
//    $location = json_decode($location);
//
//    curl_close($ch);
//
//
//
//    $loc = "";
//
//    if($location===FALSE) return "";
//
//    if (empty($location->desc)) {
//
//        $loc = $location->province.$location->city.$location->district.$location->isp;
//
//    }else{
//
//        $loc = $location->desc;
//
//    }
//
//    return $loc;

// if (0 == false) {
// 	echo 123;
// }




// $data = json_decode($data);






?>